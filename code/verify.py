"""Independent exact event enumeration, polygon/HiGHS cross-checks, deviations."""
from pathlib import Path
import sys,json,itertools,csv
import numpy as np
from model import primitives,branch,solve,central,scalar
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]

def events(p,q):
 # Independently enumerate types, signals, review/no-review and test outcomes.
 pi=p['pi'];a=p['alpha'];f=p['f'];z=p.get('zeta',1)
 types=[('G',pi,p['V']),('B1',(1-pi)*p['lam'],p['b1']),('B2',(1-pi)*(1-p['lam']),p['b2'])]
 out=dict(cost=0.,reviews=0.,replacement=0.,filing=0.,local_final=0.)
 deviations={}
 for name,mass,benefit in types:
  w=[1-a,a] if name=='G' else [a,1-a]
  submit_pay=-f
  for s in range(2):
   # Reviewed good always passes; reviewed bad passes with probability 1-z.
   submit_pay+=w[s]*benefit*(1-q[s]*(0 if name=='G' else z))
  submit=submit_pay>1e-10 # Only used away from exact ties in moment checks.
  deviations[name]=dict(submit=bool(submit),gain_submit=submit_pay,gain_decline=0.)
  if not submit:
   if name=='G':out['cost']+=mass*p['V']
   continue
  out['filing']+=mass;out['cost']+=mass*f
  for s in range(2):
   for audit,prob_a in [(0,1-q[s]),(1,q[s])]:
    for caught,prob_z in ([(0,1.)] if not audit or name=='G' else [(1,z),(0,1-z)]):
     prob=mass*w[s]*prob_a*prob_z
     action=1-caught
     out['reviews']+=prob*audit;out['replacement']+=prob*caught;out['local_final']+=prob*(1-audit)
     out['cost']+=prob*(audit*p['c']+(p['H'] if action and name!='G' else 0))
 return out,deviations

baseA=dict(pi=.5,lam=.7,H=3.,f=.2,c=.1,K=.3,b1=.3,b2=1.2,V=2.,d=.22,zeta=.9)
baseB=dict(pi=.8,lam=.5,H=3.,f=.1,c=.1,K=.1,b1=1.,b2=2.,V=1.,d=.25,zeta=.9)
rows=[];deviation_log=[]
for label,base,alphas in [('expansion',baseA,[.55,.9]),('centralization',baseB,[.55,.95])]:
 for a in alphas:
  p={**base,'alpha':a};s,branches=solve(p);slp,_=solve(p,method='lp')
  assert abs(s['cost']-slp['cost'])<1e-8
  cp=central(p);assert abs(cp-central(p,'lp'))<1e-8
  direct,dev=events(p,s['q'])
  for key in direct:assert abs(direct[key]-s[key])<1e-9,(key,direct,s)
  assert all(abs(v['gain_submit'])>1e-3 for v in dev.values())
  margin,_=solve(p,margin=.001);assert abs(margin['cost']-s['cost'])<1e-8
  D=p['f']+(1-p['pi'])*p['H'];ban=p['pi']*p['V'];choice=min({'S':s['cost'],'C':cp,'D':D,'ban':ban},key={'S':s['cost'],'C':cp,'D':D,'ban':ban}.get)
  rows.append(dict(example=label,alpha=a,selective_cost=s['cost'],central_cost=cp,delegation_cost=D,chosen=choice,**{k:s[k] for k in ['j','QB','QG','reviews','replacement','filing','local_final','reserved']}))
  deviation_log.append(dict(example=label,alpha=a,q=s['q'],types=dev))
# 160 parameter settings, 3 branches each and separate central LP.
rng=np.random.default_rng(610206);comparisons=0
for case in range(160):
 p=dict(pi=float(rng.uniform(.15,.85)),lam=float(rng.uniform(.1,.9)),H=float(rng.uniform(.5,5)),f=.1,c=float(rng.uniform(.01,.4)),K=float(rng.uniform(.01,.9)),b1=float(rng.uniform(.12,.8)),b2=float(rng.uniform(.9,2)),V=float(rng.uniform(.5,3)),d=.2,alpha=float(rng.uniform(.5,.99)),zeta=float(rng.choice([0,.6,.9,1.])))
 for j in range(3):
  x=branch(p,j);y=branch(p,j,method='lp');assert (x is None)==(y is None)
  if x is not None:
   assert abs(x['cost']-y['cost'])<1e-8
   assert abs(x['cost']-scalar(p,j))<1e-8
  comparisons+=1
 assert abs(central(p)-central(p,'lp'))<1e-8;comparisons+=1
# Dense grid is an additional upper bound, not the optimization proof.
for base in [baseA,baseB]:
 p={**base,'alpha':.9};opt,_=solve(p);best=1e9
 for qh in np.linspace(0,1,151):
  for ql in np.linspace(0,1,151):
   q=[qh,ql];_,_,m=primitives(p)
   if m@q>p['K']+1e-12:continue
   val,_=events(p,q);best=min(best,val['cost'])
 assert best>=opt['cost']-1e-8 and best-opt['cost']<.02
# Threshold ties: both pure choices and every mixture are best responses.
p={**baseA,'alpha':.9};t=(1-p['f']/p['b1'])/p['zeta'];q=[t/p['alpha'],0]
assert abs(p['b1']*(1-p['zeta']*t)-p['f'])<1e-12
# Capacity boundary: exact maximal deterrence threshold is not strict deterrence.
p0={**p,'K':.5*q[0]};assert branch(p0,1) is None
# Zero detection: no deterrence and no beneficial review.
p0={**p,'zeta':0};s,_=solve(p0);assert s['j']==0 and max(s['q'])<1e-9
# Zero capacity: no review and all submit (f < b1).
p0={**p,'K':0};s,_=solve(p0);assert abs(s['cost']-(p['f']+(1-p['pi'])*p['H']))<1e-9
# Free proposals: no strict deterrence below full review with full signal support.
p0={**p,'f':0,'zeta':1};s,_=solve(p0);assert s['j']==0
# Uninformative signal gives QB=QG under every rule.
for q in [[0,1],[.2,.7],[1,1]]:
 g,b,m=primitives({**p,'alpha':.5});assert abs(g@q-b@q)<1e-12
for folder in ['results','figures']:(ROOT/folder).mkdir(exist_ok=True)
with (ROOT/'results/examples.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
(ROOT/'results/deviations.json').write_text(json.dumps(deviation_log,indent=2))
(ROOT/'results/verification.json').write_text(json.dumps(dict(seed=610206,parameter_settings=160,optimizer_comparisons=comparisons,examples=rows,boundary_checks='passed',event_checks='passed',dense_grid_checks='passed'),indent=2))
(ROOT/'results/environment.txt').write_text('Python '+sys.version+'\nnumpy '+np.__version__+'\nscipy '+scipy.__version__+'\nmatplotlib '+matplotlib.__version__+'\n')
fig,axs=plt.subplots(1,2,figsize=(10,3.8))
for ax,base,title in zip(axs,[baseA,baseB],['A: selective initiation becomes cheaper','B: central choice becomes cheaper']):
 aa=np.linspace(.5,.98,121);ss=[];cc=[]
 for a in aa:
  p={**base,'alpha':float(a)};x,_=solve(p);ss.append(x['cost']);cc.append(central(p))
 ax.plot(aa,ss,label='Selective review infimum');ax.plot(aa,cc,label='Central choice');ax.set(xlabel='Screen accuracy alpha',ylabel='Total real cost',title=title);ax.legend(fontsize=8)
fig.suptitle('Synthetic solved examples; all costs and capacity fixed within each panel',fontsize=10);fig.tight_layout();fig.savefig(ROOT/'figures/policy_costs.png',dpi=180)
print(json.dumps(dict(comparisons=comparisons,examples=rows),indent=2))

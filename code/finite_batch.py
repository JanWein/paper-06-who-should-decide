"""R3 finite batch certificates. No downloads, random draws, or data.
Run: python3 code/finite_batch.py. Outputs results/R3_checks.json.
Policy per signal pair (HH,HL,LH,LL): u1,u2,v1,v2.
u: test probability with that sole filer; v: with both filers.
All lotteries select at most one test in every realized history.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import numpy as np
from scipy.optimize import linprog
S=list(product(range(2), repeat=2))

def weights(pi,alpha):
    B=[alpha,1-alpha]; G=[1-alpha,alpha]
    return B,G,[pi*G[k]+(1-pi)*B[k] for k in range(2)]

def coefficients(pi,alpha):
    B,G,_=weights(pi,alpha); A=[[F(0)]*16 for _ in range(2)]; T=[F(0)]*16
    for k,(s,r) in enumerate(S):
        j=4*k
        A[0][j]=(1-pi)*B[s]*B[r]; A[1][j+1]=(1-pi)*B[s]*B[r]
        A[0][j+2]=pi*B[s]*G[r]; A[1][j+3]=pi*G[s]*B[r]
        T[j]=pi*(1-pi)*G[s]*B[r]; T[j+1]=pi*(1-pi)*B[s]*G[r]
        T[j+2]=T[j+3]=pi*pi*G[s]*G[r]
    return A,T

def dot(a,b): return sum(x*y for x,y in zip(a,b))
def priority():
    x=[]
    for s,r in S:
        v=F(1,2) if s==r else F(s<r)
        x += [F(1),F(1),v,1-v]
    return x

def raw(pi,alpha,zeta,f,c,H,p,x,fixed=None):
    """Enumerate qualities, signals, filing outcomes, and test lottery.
    H is real bad-implementation harm; private benefit b is not a cost.
    Fixed assignment probabilities remain independent of filing.
    """
    B,G,_=weights(pi,alpha); out=dict(filings=F(0),tests=F(0),harm=F(0))
    for q in product(range(2),repeat=2): # 0 good, 1 bad
      probq=product_prob([1-pi if t else pi for t in q])
      for k,s in enumerate(S):
        probs=probq*product_prob([(B if q[i] else G)[s[i]] for i in range(2)])
        for fil in product(range(2),repeat=2):
          pf=probs*product_prob([(p[i] if q[i] else F(1)) if fil[i] else (1-p[i] if q[i] else F(0)) for i in range(2)])
          if not pf: continue
          if fixed is not None: audits=[fixed[k][i]*fil[i] for i in range(2)]
          elif fil==(1,1): audits=x[4*k+2:4*k+4]
          elif fil==(1,0): audits=[x[4*k],F(0)]
          elif fil==(0,1): audits=[F(0),x[4*k+1]]
          else: audits=[F(0),F(0)]
          assert sum(audits)<=1
          out['filings']+=pf*sum(fil); out['tests']+=pf*sum(audits)
          out['harm']+=pf*H*sum(fil[i]-zeta*audits[i] for i in range(2) if q[i])
    out['loss']=f*out['filings']+c*out['tests']+out['harm']
    return out

def product_prob(xs):
    z=F(1)
    for v in xs:z*=v
    return z

def interval_intersect(lo,hi,lc,hc,newlo,newhi,nlc,nhc):
    if newlo>lo:lo,lc=newlo,nlc
    elif newlo==lo:lc=lc and nlc
    if newhi<hi:hi,hc=newhi,nhc
    elif newhi==hi:hc=hc and nhc
    if lo>hi or (lo==hi and not(lc and hc)):return None
    return lo,hi,lc,hc

def all_equilibria(A,D,t):
    """Exact nine-support characterization as Cartesian intervals.
    State 0/1 pins p; state M gives open (0,1). Best-response constraint
    of each player restricts the OTHER coordinate. Handles D=0, signed D,
    equality, asymmetric equilibria and continua without grid search.
    """
    result=[]
    for support in product(['0','M','1'],repeat=2):
      box=[(F(0),F(0),True,True) if a=='0' else (F(1),F(1),True,True) if a=='1' else (F(0),F(1),False,False) for a in support]
      for i in range(2):
        # condition A_i-D_i*p_j -t >=0 for abstention, <=0 filing
        a=A[i]-t; d=D[i]; mode=support[i]
        if d==0:
          ok=(a>=0 if mode=='0' else a<=0 if mode=='1' else a==0)
          if not ok:box=None;break
          continue
        root=a/d
        if mode=='M': constraint=(root,root,True,True)
        else:
          lower=(mode=='0' and d<0) or (mode=='1' and d>0)
          constraint=(root,F(1),True,True) if lower else (F(0),root,True,True)
        j=1-i; box[j]=interval_intersect(*box[j],*constraint)
        if box[j] is None:box=None;break
      if box is not None:result.append(dict(support=support,box=box))
    return result

def target_candidates(pi,alpha,d):
    """Enumerate 12^4 hard-feasible deterministic policies. Project onto
    (mean bad-deviation scrutiny, target workload), then take lower hull.
    Label symmetrization makes the mean the common scrutiny. The optimum
    is one vertex or a mixture of two adjacent hull vertices.
    """
    A,T=coefficients(pi,alpha); av=[(a+b)/2 for a,b in zip(*A)]
    choices=[(u,v,a,b) for u,v in product([F(0),F(1)],repeat=2) for a,b in [(F(0),F(0)),(F(1),F(0)),(F(0),F(1))]]
    points={}
    for blocks in product(choices,repeat=4):
      x=sum((list(v) for v in blocks),[]); a=dot(av,x); cost=dot(T,x)
      if a not in points or cost<points[a][0]:points[a]=(cost,x)
    hull=[]
    for a,(cost,x) in sorted(points.items()):
      while len(hull)>1:
        a0,c0,_=hull[-2];a1,c1,_=hull[-1]
        if (c1-c0)*(a-a1)>=(cost-c1)*(a1-a0):hull.pop()
        else:break
      hull.append((a,cost,x))
    for left,right in zip(hull,hull[1:]):
      a,l,x=left;b,r,y=right
      if a<=d<=b:
        lam=(d-a)/(b-a);return l+lam*(r-l),hull,dict(left=a,right=b,right_weight=lam)
    raise ValueError('IC margin infeasible')

def lp_target(pi,alpha,d):
    A,T=coefficients(pi,alpha); ub=[];rhs=[]
    for k in range(4):
      row=[0]*16;row[4*k+2]=row[4*k+3]=1;ub.append(row);rhs.append(1)
    ub += [[-float(v) for v in row] for row in A];rhs += [-float(d)]*2
    z=linprog(list(map(float,T)),A_ub=ub,b_ub=rhs,bounds=[(0,1)]*16,method='highs');assert z.success
    return z.fun

def frontiers(pi,alpha):
    B,G,m=weights(pi,alpha)
    # All 81 deterministic preassignments, then LP over their mixtures.
    vals=[]
    for rule in product(range(3),repeat=4):
      q=[F(0),F(0)]
      for k,(s,r) in enumerate(S):
        if rule[k]==1:q[0]+=B[s]*m[r]
        if rule[k]==2:q[1]+=m[s]*B[r]
      vals.append(q)
    mat=[[-float(q[i]) for q in vals]+[1.] for i in range(2)]
    fit=linprog([0.]*81+[-1.],A_ub=mat,b_ub=[0,0],A_eq=[[1.]*81+[0]],b_eq=[1],bounds=[(0,None)]*82,method='highs');assert fit.success
    qf=F(1,2)+pi*(alpha-F(1,2));qr=1-pi+pi*alpha
    assert abs(fit.x[-1]-float(qf))<1e-9
    # Independent raw-probability recourse LP, 16 policy coordinates plus q.
    A,_=coefficients(pi,alpha);ub=[[-float(v) for v in a]+[1] for a in A];rhs=[0,0]
    for k in range(4):
      row=[0]*17;row[4*k+2]=row[4*k+3]=1;ub.append(row);rhs.append(1)
    fit=linprog([0.]*16+[-1.],A_ub=ub,b_ub=rhs,bounds=[(0,1)]*16+[(0,None)],method='highs');assert fit.success
    assert abs(fit.x[-1]-float(qr))<1e-9
    return dict(pi=str(pi),alpha=str(alpha),fixed=str(qf),recourse=str(qr))

def main():
    checks={}; checks['frontier_cases']=[frontiers(p,a) for p in [F(1,5),F(1,2),F(4,5),F(1)] for a in [F(1,2),F(3,5),F(7,10),F(4,5),F(9,10),F(1)]]
    pi=alpha=F(1,2);z=F(9,10);f=F(2,5);c=F(1,10);H=F(3);t=F(2,3);x=priority()
    checks['witness']={}
    for p,loss in [(F(0),F(19,40)),(F(1,3),F(46,45)),(F(1),F(51,20))]:
      r=raw(pi,alpha,z,f,c,H,[p,p],x);assert r['loss']==loss;checks['witness'][str(p)]=r
    fixed=[(F(7,10),F(3,10))]*4
    r=raw(pi,alpha,z,f,c,H,[F(0),F(1)],x,fixed);assert r['loss']==F(44,25);checks['fixed_witness']=r
    assert 1-z*F(3,4)-f==-F(3,40)
    # Boundary cases and asymmetries through exact support solver.
    checks['priority_equilibria']={}
    for t0,count in [(F(2,5),1),(F(1,2),2),(F(2,3),3),(F(3,4),2),(F(4,5),1)]:
      eq=all_equilibria([F(3,4)]*2,[F(1,4)]*2,t0);assert len(eq)==count;checks['priority_equilibria'][str(t0)]=eq
    assert all_equilibria([F(7,10),F(3,10)],[F(0)]*2,t)[0]['support']==('0','1')
    assert len(all_equilibria([F(1,4)]*2,[-F(1,2)]*2,F(1,2)))==3 # negative D: two asymmetric pure BNE and interior mix
    assert len(all_equilibria([t,t],[F(0)]*2,t))==9 # entire square partitioned into nine supports
    checks['target_cost']={}
    for a,d in [(F(1,2),F(41,60)),(F(3,4),F(41,60))]:
      val,hull,mix=target_candidates(pi,a,d); independent=lp_target(pi,a,d);assert abs(float(val)-independent)<1e-9
      checks['target_cost'][str(a)]=dict(d=d,tests=val,hull=[(a,c) for a,c,_ in hull],mixture=mix)
    assert checks['target_cost']['1/2']['tests']==F(41,60)
    assert f+c*F(41,60)==F(281,600)
    # Exact affine certificate for the informative region: at upper excluded
    # endpoint bound=3/5, slope=-27/20, hence strictly larger below endpoint.
    assert F(69,40)-F(27,20)*F(5,6)==F(3,5)>F(19,40)
    assert F(17,20)+F(13,10)*F(2,3)==F(103,60)
    assert F(103,60)+F(13,10)*F(1,30)==F(44,25)
    checks['informative_region']=dict(interval='[1/2,5/6)',harm_lower_bound='69/40 - 27 alpha/20',endpoint='3/5',target_loss='19/40')
    checks['status']='all deterministic assertions passed; no empirical data or random draws'
    out=Path(__file__).resolve().parents[1]/'results'/'R3_checks.json';out.write_text(json.dumps(checks,default=str,indent=2)+'\n');print(out);print(checks['status'])
if __name__=='__main__':main()

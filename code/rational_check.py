"""Exact rational certificates for both main examples; standard library only."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def vertices(A,z):
 out=[]
 for i,j in combinations(range(len(z)),2):
  a,b=A[i];c,d=A[j];det=a*d-b*c
  if det==0:continue
  x=((z[i]*d-b*z[j])/det,(a*z[j]-z[i]*c)/det)
  if all(dot(row,x)<=bound for row,bound in zip(A,z)):out.append(x)
 return out

def certificate(p):
 pi,l,H,V,f,c,d,K,b1,b2,a,z=p
 g=(1-a,a);b=(a,1-a);m=tuple(pi*x+(1-pi)*y for x,y in zip(g,b));t=[F(0),(1-f/b1)/z,(1-f/b2)/z,F(1)]
 square=[(-F(1),F(0)),(F(0),-F(1)),(F(1),F(0)),(F(0),F(1)),m];bounds=[F(0),F(0),F(1),F(1),K]
 Qmax=max(dot(b,q) for q in vertices(square,bounds))
 branches=[]
 for j,B in enumerate([1-pi,(1-pi)*(1-l),F(0)]):
  if j and Qmax<=t[j]:continue
  A=square+[b,tuple(-x for x in b)];rhs=bounds+[min(t[j+1],F(1)),-t[j]]
  qs=vertices(A,rhs)
  if not qs:continue
  def loss(q):return f*(pi+B)+c*pi*dot(g,q)+(c-H*z)*B*dot(b,q)+H*B
  q=min(qs,key=loss);Q=dot(b,q)
  branches.append((loss(q),j,q,Q,B))
 L,j,q,Q,B=min(branches,key=lambda x:x[0]);assert j==0 or Q>t[j]
 gains=[V-f,b1*(1-z*Q)-f,b2*(1-z*Q)-f]
 assert all(x!=0 for x in gains)
 # Central geometric problem with exact negative-test updating.
 e=tuple(min(pi*x*V,(1-pi)*y*H) for x,y in zip(g,b))
 aud=tuple(c*ms+min(pi*x*V,(1-z)*(1-pi)*y*H) for ms,x,y in zip(m,g,b))
 vals=[sum(e)+dot(tuple(x-y for x,y in zip(aud,e)),r)+d for r in vertices(square,bounds)]
 C=min(vals);D0=pi*V;D1=f+(1-pi)*H
 costs={'S':L,'C':C,'D0':D0,'D1':D1};choice=min(costs,key=costs.get)
 return dict(q=[str(x) for x in q],QB=str(Q),S_cost=str(L),C_cost=str(C),D0_cost=str(D0),D1_cost=str(D1),choice=choice,branch=j,gains=[str(x) for x in gains],replacement=str(z*B*Q),all_branch_costs=[str(x[0]) for x in branches])
cases={}
for label,base,alphas in [('expansion',['.5','.7','3','2','.2','.1','.22','.3','.3','1.2'],['.55','.9']),('centralization',['.8','.5','3','1','.1','.1','.25','.1','1','2'],['.55','.95'])]:
 for a in alphas:
  cases[label+'_'+a]=certificate(tuple(F(x) for x in base+[a,'.9']))
assert [v['choice'] for v in cases.values()]==['C','S','S','C']
(ROOT/'results/rational_certificates.json').write_text(json.dumps(cases,indent=2))
print('Four exact rational optimizer certificates passed.')

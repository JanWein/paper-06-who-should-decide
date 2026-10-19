"""Finite-type strategic screening, exact polygon enumeration and independent LP."""
import itertools
import numpy as np
from scipy.optimize import linprog

def primitives(p):
 pi=p['pi'];a=p['alpha'];g=np.array([1-a,a]);b=np.array([a,1-a]);m=pi*g+(1-pi)*b
 return g,b,m

def branch(p,j,margin=0,method='vertices'):
 g,b,m=primitives(p);pi=p['pi'];f=p['f'];c=p['c'];H=p['H'];K=p['K']
 zeta=p.get('zeta',1.0)
 # Robust branch j>0 exists only if strict deterrence is feasible.
 qmaxH=min(1,K/m[0]);qmaxL=min(1,max(0,K-m[0]*qmaxH)/m[1]);maxQB=float(b@np.array([qmaxH,qmaxL]))
 thresholds=[0,(1-f/p['b1'])/zeta if zeta else 2,(1-f/p['b2'])/zeta if zeta else 3,1]
 if j and maxQB<=thresholds[j]+1e-10:return None
 B=[1-pi,(1-pi)*(1-p['lam']),0][j]
 lo=thresholds[j]+(margin if j else 0);hi=min(1,thresholds[j+1]-(margin if j<2 and thresholds[j+1]<=1 else 0))
 if lo>hi:return None
 A=np.array([[-1,0],[0,-1],[1,0],[0,1],m,b,-b]);z=np.array([0,0,1,1,K,hi,-lo])
 obj=c*pi*g+(c-H*zeta)*B*b;const=f*(pi+B)+H*B
 if method=='lp':
  res=linprog(obj,A_ub=A,b_ub=z,bounds=[(None,None)]*2,method='highs')
  if not res.success:return None
  q=res.x
 else:
  candidates=[]
  for i,k in itertools.combinations(range(len(z)),2):
   M=A[[i,k]]
   if abs(np.linalg.det(M))<1e-12:continue
   x=np.linalg.solve(M,z[[i,k]])
   if np.all(A@x<=z+1e-9):candidates.append(x)
  if not candidates:return None
  q=min(candidates,key=lambda x:obj@x)
 QG=float(g@q);QB=float(b@q);reviews=pi*QG+B*QB
 return dict(cost=float(const+obj@q),j=j,q=q.tolist(),QB=QB,QG=QG,bad_mass=B,replacement=zeta*B*QB,reviews=reviews,reserved=float(m@q),local_final=pi+B-reviews,filing=pi+B)

def solve(p,margin=0,method='vertices'):
 rows=[branch(p,j,margin,method) for j in range(3)];rows=[x for x in rows if x is not None]
 return min(rows,key=lambda x:x['cost']),rows

def central(p,method='vertices'):
 g,b,m=primitives(p);pi=p['pi'];V=p['V'];H=p['H'];c=p['c'];K=p['K'];d=p.get('d',0)
 # per signal choose accept/reject/review probabilities; independent LP benchmark
 zeta=p.get('zeta',1.0)
 costs=np.array([(1-pi)*b*H,pi*g*V,c*m+np.minimum(pi*g*V,(1-pi)*b*(1-zeta)*H)]).T
 if method=='lp':
  obj=costs.ravel();Aeq=np.array([[1,1,1,0,0,0],[0,0,0,1,1,1]])
  cap=np.array([0,0,m[0],0,0,m[1]])
  res=linprog(obj,A_ub=[cap],b_ub=[K],A_eq=Aeq,b_eq=[1,1],bounds=(0,1),method='highs')
  assert res.success
  return float(res.fun+d)
 # C value = min over review amounts r of error*(1-r)+c*m*r.
 err=np.minimum(costs[:,0],costs[:,1]);coef=costs[:,2]-err
 A=np.array([[-1,0],[0,-1],[1,0],[0,1],m]);z=np.array([0,0,1,1,K]);vs=[]
 for i,k in itertools.combinations(range(5),2):
  if abs(np.linalg.det(A[[i,k]]))<1e-12:continue
  r=np.linalg.solve(A[[i,k]],z[[i,k]])
  if np.all(A@r<=z+1e-9):vs.append(float(err.sum()+coef@r+d))
 return min(vs)


def scalar(p,j):
 """Independent one-dimensional frontier solver; closed branch infimum."""
 a=p['alpha'];pi=p['pi'];f=p['f'];z=p.get('zeta',1);K=p['K']
 mh=pi*(1-a)+(1-pi)*a
 qm=a*K/mh if K<=mh else a+(1-a)*(K-mh)/(1-mh)
 qm=min(1,qm)
 t=[0,(1-f/p['b1'])/z if z else 2,(1-f/p['b2'])/z if z else 3,1]
 if j and qm<=t[j]+1e-10:return None
 lo=t[j];hi=min(1,t[j+1],qm)
 if lo>hi:return None
 B=[1-pi,(1-pi)*(1-p['lam']),0][j]
 gamma=lambda Q:(1-a)*Q/a if Q<=a else 1-a+a*(Q-a)/(1-a)
 loss=lambda Q: f*(pi+B)+p['c']*pi*gamma(Q)+(p['c']-p['H']*z)*B*Q+p['H']*B
 qs=[lo,hi]+([a] if lo<=a<=hi else [])
 Q=min(qs,key=loss)
 return float(loss(Q))

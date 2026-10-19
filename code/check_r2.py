"""Focused R2 regression: exact fractions, deterministic inputs, no random sweep."""
from fractions import Fraction as F
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
def gamma(a,Q):
 return (1-a)*Q/a if Q<=a else 1-a+a*(Q-a)/(1-a)
rows=[]
for Q in [F(3,5),F(7,10),F(9,10)]:
 assert gamma(Q,Q)==1-Q
 left_limit=-1/(1-Q);right_limit=-1/Q
 assert left_limit!=right_limit and left_limit<0 and right_limit<0
 previous=None
 for h in [F(1,100),F(1,10000),F(1,1000000)]:
  left=(gamma(Q,Q)-gamma(Q-h,Q))/h
  right=(gamma(Q+h,Q)-gamma(Q,Q))/h
  errors=(abs(left-left_limit),abs(right-right_limit))
  if previous:assert errors[0]<previous[0] and errors[1]<previous[1]
  previous=errors
 assert errors[0]<F(1,1000) and errors[1]<F(1,1000)
 rows.append(dict(Q=str(Q),value=str(gamma(Q,Q)),left_limit=str(left_limit),right_limit=str(right_limit),left_difference=float(left),right_difference=float(right)))
# Exact monotonicity checks on a deterministic mesh across and away from the kink.
for Q in [F(1,5),F(1,2),F(7,10),F(9,10)]:
 vals=[gamma(F(j,100),Q) for j in range(50,100)]
 assert all(x>y for x,y in zip(vals,vals[1:]))
assert gamma(F(1,2),F(1,2))==F(1,2)
assert -1/(1-F(1,2))==-1/F(1,2)==-2
for a in [F(1,2),F(7,10),F(9,10)]:
 assert gamma(a,F(0))==0 and gamma(a,F(1))==1
pi=F(1,2);B=F(3,20);QB=F(27,50);QG=F(3,50)
R=pi*QG+(1-pi)*QB;J=pi*QG+B*QB
assert R==F(3,10) and J==F(111,1000)
assert R-J==((1-pi)-B)*QB==F(189,1000)
result=dict(status='passed',one_sided=rows,reservation=str(R),workload=str(J),gap=str(R-J),random_draws=0)
(ROOT/'results/R2_regression.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))

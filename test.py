from alphas import *

Qcd(0)
t=0.42
print(  "RG-solver gives : "+str( Solver(t,2) )  )
print(  "Exact solution  : "+str( A_2loop(t ) )  )

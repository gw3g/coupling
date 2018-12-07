from alphas import *

out = open("error2.dat",'w')
out.write("# Columns: t, UV, a1(tinf=1e4), a2(tinf=1e2)\n")
t_s,a_U=solver(.5,1e2,50,1e4)
t_s,a_L=solver(.5,1e2,50,1e3)
for i in range(len(t_s)):
    out.write( "{0:.5e}  {1:.5e}  {2:.5e}  {3:.5e}\n".format(t_s[i],a2loop(t_s[i]),a_U[i],a_L[i]) )
out.close()

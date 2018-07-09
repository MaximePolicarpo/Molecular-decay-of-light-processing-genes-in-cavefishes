#!/usr/bin/env python

import math
import matplotlib.pyplot as plt

# D = Number of pseudoenes
# T = Number of total genes
# U = Lof frequency(lof) * mean gene length(l) * mutation rate(u)
# t = time in generations

#T=85 #Astyanax
T=76 #Lucifuga

l = 1080.0 #Both Astyanax and Lucifuga
u = 1e-8
lof= 0.086 #Stop = 0.04, FS = 0.04, Splice = 0.02, LossStartorStop = 0.006


t_fs=list()
for D in range(0,T):
	t_fs.append((-1/(u*lof*l))*math.log((float(T)-float(D))/float(T)))


lof = 0.04
t_s=list()
for D in range(0,T):
	t_s.append((-1/(u*lof*l))*math.log((float(T)-float(D))/float(T)))

lof = 0.106
t_fss=list()
for D in range(0,T):
	t_fss.append((-1/(u*lof*l))*math.log((float(T)-float(D))/float(T)))



print(t_fss[19])


y=range(0, T)
plt.plot(t_s, y, label="Stop (0.04*10-8)")
plt.plot(t_fs, y, label="Stop+Frameshift (0.086*10-8)")
plt.plot(t_fss, y, label="Stop+Frameshift+Splice (0.106*10-8)")
plt.xlabel("Number of generations")
plt.ylabel("Number of fixed pseudogenes")
plt.legend()




plt.show()





#print((-1/(u*lof*l))*math.log((85.0-1.0)/85.0))

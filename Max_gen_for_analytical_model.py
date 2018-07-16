#!/usr/bin/env python

#Corresponding auhors : Didier Casane (didier.casane@egce.cnrs-gif.fr) and Maxime Policarpo (maxime.policarpo@hotmail.fr)

import math
import matplotlib.pyplot as plt

T=85 #number of genes for Astyanax mexicanus
#T=76 #number of genes for Lucifuga dentata

l=1061 #mean size of Astyanax mexicanus genes
#l = 1091 #mean size of Lucifuga dentata genes
u = 1e-8
lof = 0.112

T_Lof=list()
for D in range(0,T):
	T_Lof.append((-1/(u*lof*l))*math.log((float(T)-float(D))/float(T)))


print(T_Lof[19])
print(T_Lof[1])


#y=range(0, T)
#plt.plot(T_Lof, y, label="u=1e-08")
#plt.plot(T2_Lof, y, label="u=3e-09")
#plt.xlim(0,50000)
#plt.ylim(0,5)
#plt.xlabel("Number of generations")
#plt.ylabel("Number of fixed pseudogenes")
#plt.legend()
#plt.show()

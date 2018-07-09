#!/usr/bin/env python3


import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


T=76   #Nombre de genes de lucifuga(76)
l=1091 #Taille moyenne des genes de lucifuga
#T=85    #Nombre de genes de Astyanax cavefish(85)
#l=1061 #taille moyenne genes astyanax

dimension = (100000, (T+1))
bino_matrice = np.zeros(dimension)



for D in range(0,(T+1)):
	myvector=list()
	for time in range(1,1000001):
        	if(time % 10 == 0):
			X1=(((math.factorial(T))/((math.factorial(D)*(math.factorial(T-D))))))*((1-math.exp((-10**-8)*0.112*time*l))**D)*(math.exp((-10**-8)*0.112*time*l))**(T-D) #FS+Stop+Introns+Loss
#			X1=(((math.factorial(T))/((math.factorial(D)*(math.factorial(T-D))))))*((1-math.exp((-7*(10**-9))*0.112*time*l))**D)*(math.exp((-7*(10**-9))*0.112*time*l))**(T-D) #FS+Stop+Introns+Loss
			myvector.append(X1)
	bino_matrice[:,D]=myvector



#print(myvector.index(np.amax(myvector)))
#print(np.amax(myvector))


print(bino_matrice)

np.savetxt(("{}genes_0.112_7e9u_1080CDS.csv".format(T)), bino_matrice, delimiter=",")




#plt.plot(myvector)
#plt.xlim(0,100000)
#plt.ylim(0,1)
#plt.show()






#with PdfPages("test_binomial_law.pdf") as pdf:
#	for i in range(77):
#		plt.plot(bino_matrice[:,i], 'k-')
#		plt.ylim([-0.05,1.05])
#       	plt.title("Probability of {} pseudogenes".format(i))
#		pdf.savefig()
#





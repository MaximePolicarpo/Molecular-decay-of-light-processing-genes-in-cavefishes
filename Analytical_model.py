#!/usr/bin/env python3

#Corresponding auhors : Didier Casane (didier.casane@egce.cnrs-gif.fr) and Maxime Policarpo (maxime.policarpo@hotmail.fr)


import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


T=76   #Number of Lucifuga dentata genes
l=1091 #mean size of Lucifuga dentata genes
#T=85    #Number of Astyanax mexicanus genes
#l=1061 #mean size of Astyanax mexicanus cavefish genes


dimension = (100000, (T+1))
bino_matrice = np.zeros(dimension)





for D in range(0,(T+1)):
	myvector=list()
	for time in range(1,1000001):
        	if(time % 10 == 0):
			X1=(((math.factorial(T))/((math.factorial(D)*(math.factorial(T-D))))))*((1-math.exp((-10**-8)*0.112*time*l))**D)*(math.exp((-10**-8)*0.112*time*l))**(T-D)
			myvector.append(X1)
	bino_matrice[:,D]=myvector



#print(myvector.index(np.amax(myvector)))
#print(np.amax(myvector))


print(bino_matrice)

np.savetxt(("Matrice_loi_binomiale.csv"), bino_matrice, delimiter=",")




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









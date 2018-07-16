###  Simulation of mutation distribution observed vs simulated #######

#Corresponding auhors : Didier Casane (didier.casane@egce.cnrs-gif.fr) and Maxime Policarpo (maxime.policarpo@hotmail.fr)
#library("psych", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.0")

# We observe 19 pseudogenes in dentata for a total of 24 mutations
V1=c(57,15,3,1)
Probas_obs=c(57/76,15/76,3/76,1/76)

#Length of our 76 genes :
#genes=c(1062,1163,1749,1672,1305,1083,1044,1200,1350,1215,1201,1063,861,1218,1014,1011,1049,888,882,1059,1038,1012,1058,1134,1143,912,957,1227,1206,1164,527,600,592,591,696,660,663,630,624,552,1596,1605,1064,1197,558,571,594,573,555,570,579,1689,1663,1665,1635,2569,2565,2568,267,264,219,261,612,582,580,3212,3321,3315,1054,1053,1024,1023,1022,810,222,213)

genes_length=c(1062,1164,1749,1672,1305,1083,1044,1200,1350,1215,1200,1062,861,1218,1014,1011,1049,888,882,1059,1038,1011,1059,1134,1143,912,957,1227,1206,1164,527,600,591,591,696,660,663,630,624,552,1596,1605,1062,1197,558,570,594,570,555,570,579,1689,1663,1665,1635,2568,2565,2568,267,264,219,261,612,582,579,3212,3321,3315,1053,1053,1023,1023,1023,810,222,213)
genes_intron=c(4,3,10,11,8,7,4,5,4,4,4,3,3,3,4,3,3,7,6,0,4,6,4,4,3,3,3,3,2,4,2,4,3,4,4,4,4,4,5,4,13,13,13,14,6,3,3,3,3,3,3,6,6,3,3,22,22,21,2,2,1,1,2,2,2,21,16,17,7,7,8,8,7,6,1,1)

genes=c()
for(i in 1:76){
  g=paste0("X", i)
  genes=c(genes, g)
}

madf=data.frame(genes, genes_length, genes_intron)

Graphique_proba=c()
Graphique_proba=rbind(Graphique_proba, Probas_obs)

Graphique_nb=c()
Graphique_nb=rbind(Graphique_nb, V1)


for(x in c(76,70,60,50,40,30,20,19)){               #We want to simulate what would happen when we reduce the number of neutral genes
  Total_Mutations=c()
  X=x
  for(i in 1:10000){
    set_de_genes=sample(madf$genes, X, replace=FALSE)  #We randomly pick x genes among our total set of genes
    probabilitees=c()
    mes_probas=c()
    proba_introns=c()
    mes_probas_introns=c()
    for(i in set_de_genes){
      probabilitees=c(probabilitees,madf[which(madf$genes == i),]$genes_length)
      proba_introns=c(proba_introns,madf[which(madf$genes == i),]$genes_intron)
    }
    for(i in probabilitees){
      mes_probas=c(mes_probas, i/sum(probabilitees))
    }
    for(i in proba_introns){
      mes_probas_introns=c(mes_probas_introns, i/sum(proba_introns))
    }
    
    simulation=c(sample(set_de_genes, 18, prob=mes_probas, replace=TRUE), sample(set_de_genes, 6, prob=mes_probas_introns, replace=TRUE))    #We randomly place 18 mutations (Frameshifts and Stops) and 6 intronic mutations
    
    mes_mutations=as.vector(table(simulation))
    deduction_0=76-length(mes_mutations)
    nb_0=rep(0,deduction_0)
    mes_mutations=c(mes_mutations, nb_0)                                   
    Total_Mutations=c(Total_Mutations, mes_mutations)
  }
  Proba_0mutation=length(Total_Mutations[Total_Mutations==0])/length(Total_Mutations)     
  Proba_1mutation=length(Total_Mutations[Total_Mutations==1])/length(Total_Mutations)     
  Proba_2mutation=length(Total_Mutations[Total_Mutations==2])/length(Total_Mutations)     
  Proba_3mutation=length(Total_Mutations[Total_Mutations==3])/length(Total_Mutations)     
  Proba_simu=c(Proba_0mutation,Proba_1mutation,Proba_2mutation,Proba_3mutation)           
  Graphique_proba=rbind(Graphique_proba, Proba_simu)  
  
  nb_0mutation=(length(Total_Mutations[Total_Mutations==0])/length(Total_Mutations))*76    
  nb_1mutation=(length(Total_Mutations[Total_Mutations==1])/length(Total_Mutations))*76     
  nb_2mutation=(length(Total_Mutations[Total_Mutations==2])/length(Total_Mutations))*76     
  nb_3mutation=(length(Total_Mutations[Total_Mutations==3])/length(Total_Mutations))*76     
  nb_simu=c(nb_0mutation,nb_1mutation,nb_2mutation,nb_3mutation)         
  Graphique_nb=rbind(Graphique_nb, nb_simu)
  
}

#Barplots 

#Probability calculated with a simple binomial law (= simulation with every genes having the same length)
proba_de_pas_muter=c()
proba_de_muter=c()
for(x in c(76,70,60,50,40,30,20,19)){  
  proba_de_pas_muter=c(proba_de_pas_muter, (1-(1/x)))
  proba_de_muter=c(proba_de_muter, (1/x))
}


proba_0_muta_per_gene=c()
proba_1_muta_per_gene=c()
proba_2_muta_per_gene=c()
proba_3_muta_per_gene=c()
nb=c(76,70,60,50,40,30,20,19)
matrice=cbind(nb, proba_de_muter)

for(x in 1:8){
  proba_0_muta_per_gene=c(proba_0_muta_per_gene, (dbinom(0,24,matrice[x,2])*(matrice[x,1]/76))+((76-matrice[x,1]))/76)
  proba_1_muta_per_gene=c(proba_1_muta_per_gene, (dbinom(1,24,matrice[x,2])*(matrice[x,1]/76)))
  proba_2_muta_per_gene=c(proba_2_muta_per_gene, (dbinom(2,24,matrice[x,2])*(matrice[x,1]/76)))
  proba_3_muta_per_gene=c(proba_3_muta_per_gene, (dbinom(3,24,matrice[x,2])*(matrice[x,1]/76)))
}



nb_0=proba_0_muta_per_gene*76
nb_1=proba_1_muta_per_gene*76
nb_2=proba_2_muta_per_gene*76
nb_3=proba_3_muta_per_gene*76

binom_distrib=cbind(proba_0_muta_per_gene, proba_1_muta_per_gene, proba_2_muta_per_gene, proba_3_muta_per_gene)
binom_distrib_nb=cbind(nb_0,nb_1,nb_2,nb_3)

binom_distrib=rbind(Probas_obs, binom_distrib)
binom_distrib_nb=rbind(V1, binom_distrib_nb)


#### Version two : Observation plotted as a line insteal of a bar #######

Graphique_nb_v2=Graphique_nb[-1,]

pdf("Figure6.pdf", width = 12, height =9)

xx <- barplot(Graphique_nb_v2, beside=T,ylim=c(0,76),col=c("gray67", "goldenrod4", "deepskyblue4", "darkorchid1", "darkolivegreen2", "coral", "cadetblue2", "bisque1"), names.arg=c(0,1,2,3))
text(x= xx,y=Graphique_nb_v2+1, label= round(Graphique_nb_v2, digits=2), pos=3, cex=0.7, col=c("gray67", "goldenrod4", "deepskyblue4", "darkorchid1", "darkolivegreen4", "coral", "cadetblue4", "bisque3"))
text(x=31,y=68, "Number of", cex=2)
text(x=31, y=64, "neutral genes", cex=2)
legend(x=29, y=62,legend=c("76", "70","60","50","40","30","20","16"), col=c("gray67", "goldenrod4", "deepskyblue4", "darkorchid1", "darkolivegreen2", "coral", "cadetblue2", "bisque1"), pch=15, cex=1.5)

x=1.5
for(c in nb_0){
  points(x=x, y=c, pch=19, col="black")
  x=x+1
}
x=10.5
for(c in nb_1){
  points(x=x, y=c, pch=19, col="black")
  x=x+1
}
x=19.5
for(c in nb_2){
  points(x=x, y=c, pch=19, col="black")
  x=x+1
}
x=28.5
for(c in nb_3){
  points(x=x, y=c, pch=19, col="black")
  x=x+1
}

lines(x=c(1,9), y=c(57,57), col="red", lty=1, lwd=3)
lines(x=c(10,18), y=c(15,15), col="red", lty=1, lwd=3)
lines(x=c(19,27), y=c(3,3), col="red", lty=1, lwd=3)
lines(x=c(28,36), y=c(1,1), col="red", lty=1, lwd=3)
text(x=9.7, y=57, col="red", "57", cex=1.5)
text(x=18.7, y=15, col="red", "15", cex=1.5)
text(x=27.7, y=3, col="red", "3", cex=1.5)
text(x=36.7, y=1, col="red", "1", cex=1.5)

title(main="Distribution of the number of LoF mutations per gene", xlab="Number of LoF mutations per gene", cex.lab=1.8, cex.main=2)
mtext("Number of genes", side=2, line=2, cex=2)
#legend(x=20, y=75, legend="Observed values", col="red", lty=1, bty="n", lwd=3, cex=1.5)

dev.off()



#########################################################################################################################################
#########################################################################################################################################
########################################### Same script but for holguinensis datas ######################################################
#########################################################################################################################################
#########################################################################################################################################


###  Simulation of mutation distribution observed vs simulated #######
library("psych", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.0")
# We observe 19 pseudogenes in dentata for a total of 24 mutations
V1=c(70,4,2,0)
Probas_obs=c(70/76,4/76,2/76,0/76)



#Length of our 76 genes :
#genes=c(1062,1163,1749,1672,1305,1083,1044,1200,1350,1215,1201,1063,861,1218,1014,1011,1049,888,882,1059,1038,1012,1058,1134,1143,912,957,1227,1206,1164,527,600,592,591,696,660,663,630,624,552,1596,1605,1064,1197,558,571,594,573,555,570,579,1689,1663,1665,1635,2569,2565,2568,267,264,219,261,612,582,580,3212,3321,3315,1054,1053,1024,1023,1022,810,222,213)

genes_length=c(1062,1164,1749,1672,1305,1083,1044,1200,1350,1215,1200,1062,861,1218,1014,1011,1049,888,882,1059,1038,1011,1059,1134,1143,912,957,1227,1206,1164,527,600,591,591,696,660,663,630,624,552,1596,1605,1062,1197,558,570,594,570,555,570,579,1689,1663,1665,1635,2568,2565,2568,267,264,219,261,612,582,579,3212,3321,3315,1053,1053,1023,1023,1023,810,222,213)
genes_intron=c(4,3,10,11,8,7,4,5,4,4,4,3,3,3,4,3,3,7,6,0,4,6,4,4,3,3,3,3,2,4,2,4,3,4,4,4,4,4,5,4,13,13,13,14,6,3,3,3,3,3,3,6,6,3,3,22,22,21,2,2,1,1,2,2,2,21,16,17,7,7,8,8,7,6,1,1)

genes=c()
for(i in 1:76){
  g=paste0("X", i)
  genes=c(genes, g)
}

madf=data.frame(genes, genes_length, genes_intron)

Graphique_proba=c()
Graphique_proba=rbind(Graphique_proba, Probas_obs)

Graphique_nb=c()
Graphique_nb=rbind(Graphique_nb, V1)


for(x in c(76,65,55,45,35,25,15,6)){               #We want to simulate what would happen when we reduce the number of neutral genes
  Total_Mutations=c()
  X=x
  for(i in 1:10000){
    set_de_genes=sample(madf$genes, X, replace=FALSE)  #We randomly pick x genes among our total set of genes
    probabilitees=c()
    mes_probas=c()
    proba_introns=c()
    mes_probas_introns=c()
    for(i in set_de_genes){
      probabilitees=c(probabilitees,madf[which(madf$genes == i),]$genes_length)
      proba_introns=c(proba_introns,madf[which(madf$genes == i),]$genes_intron)
    }
    for(i in probabilitees){
      mes_probas=c(mes_probas, i/sum(probabilitees))
    }
    for(i in proba_introns){
      mes_probas_introns=c(mes_probas_introns, i/sum(proba_introns))
    }
    
    simulation=c(sample(set_de_genes, 3, prob=mes_probas, replace=TRUE), sample(set_de_genes, 5, prob=mes_probas_introns, replace=TRUE))    #We randomly place 18 mutations (Frameshifts and Stops) and 6 intronic mutations
    
    mes_mutations=as.vector(table(simulation))
    deduction_0=76-length(mes_mutations)
    nb_0=rep(0,deduction_0)
    mes_mutations=c(mes_mutations, nb_0)                                   
    Total_Mutations=c(Total_Mutations, mes_mutations)
  }
  Proba_0mutation=length(Total_Mutations[Total_Mutations==0])/length(Total_Mutations)     
  Proba_1mutation=length(Total_Mutations[Total_Mutations==1])/length(Total_Mutations)     
  Proba_2mutation=length(Total_Mutations[Total_Mutations==2])/length(Total_Mutations)     
  Proba_3mutation=length(Total_Mutations[Total_Mutations==3])/length(Total_Mutations)     
  Proba_simu=c(Proba_0mutation,Proba_1mutation,Proba_2mutation,Proba_3mutation)           
  Graphique_proba=rbind(Graphique_proba, Proba_simu)  
  
  nb_0mutation=(length(Total_Mutations[Total_Mutations==0])/length(Total_Mutations))*76    
  nb_1mutation=(length(Total_Mutations[Total_Mutations==1])/length(Total_Mutations))*76     
  nb_2mutation=(length(Total_Mutations[Total_Mutations==2])/length(Total_Mutations))*76     
  nb_3mutation=(length(Total_Mutations[Total_Mutations==3])/length(Total_Mutations))*76     
  nb_simu=c(nb_0mutation,nb_1mutation,nb_2mutation,nb_3mutation)         
  Graphique_nb=rbind(Graphique_nb, nb_simu)
  
}

#Barplots 

#Probability calculated with a simple binomial law (= simulation with every genes having the same length)
proba_de_pas_muter=c()
proba_de_muter=c()
for(x in c(76,65,55,45,35,25,15,6)){  
  proba_de_pas_muter=c(proba_de_pas_muter, (1-(1/x)))
  proba_de_muter=c(proba_de_muter, (1/x))
}


proba_0_muta_per_gene=c()
proba_1_muta_per_gene=c()
proba_2_muta_per_gene=c()
proba_3_muta_per_gene=c()
nb=c(76,65,55,45,35,25,15,6)
matrice=cbind(nb, proba_de_muter)

for(x in 1:8){
  proba_0_muta_per_gene=c(proba_0_muta_per_gene, (dbinom(0,8,matrice[x,2])*(matrice[x,1]/76))+((76-matrice[x,1]))/76)
  proba_1_muta_per_gene=c(proba_1_muta_per_gene, (dbinom(1,8,matrice[x,2])*(matrice[x,1]/76)))
  proba_2_muta_per_gene=c(proba_2_muta_per_gene, (dbinom(2,8,matrice[x,2])*(matrice[x,1]/76)))
  proba_3_muta_per_gene=c(proba_3_muta_per_gene, (dbinom(3,8,matrice[x,2])*(matrice[x,1]/76)))
}



nb_0=proba_0_muta_per_gene*76
nb_1=proba_1_muta_per_gene*76
nb_2=proba_2_muta_per_gene*76
nb_3=proba_3_muta_per_gene*76

binom_distrib=cbind(proba_0_muta_per_gene, proba_1_muta_per_gene, proba_2_muta_per_gene, proba_3_muta_per_gene)
binom_distrib_nb=cbind(nb_0,nb_1,nb_2,nb_3)

binom_distrib=rbind(Probas_obs, binom_distrib)
binom_distrib_nb=rbind(V1, binom_distrib_nb)


#### Version two : Observation plotted as a line insteal of a bar #######

Graphique_nb_v2=Graphique_nb[-1,]

pdf("Figure6.pdf", width = 12, height =9)

xx <- barplot(Graphique_nb_v2, beside=T,ylim=c(0,76),col=c("gray67", "goldenrod4", "deepskyblue4", "darkorchid1", "darkolivegreen2", "coral", "cadetblue2", "bisque1"), names.arg=c(0,1,2,3), cex.names=1.5, cex.axis = 1.5)
text(x= xx,y=Graphique_nb_v2+1, label= round(Graphique_nb_v2, digits=0), pos=3, cex=1.3, col=c("gray67", "goldenrod4", "deepskyblue4", "darkorchid1", "darkolivegreen4", "coral", "cadetblue4", "bisque3"))
text(x=31,y=68, "Number of", cex=2)
text(x=31, y=64, "neutral genes", cex=2)
legend(x=29, y=62,legend=c("76","65","55","45","35","25","15","6"), col=c("gray67", "goldenrod4", "deepskyblue4", "darkorchid1", "darkolivegreen2", "coral", "cadetblue2", "bisque1"), pch=15, cex=1.5)

x=1.5
for(c in nb_0){
  points(x=x, y=c, pch=19, col="black")
  x=x+1
}
x=10.5
for(c in nb_1){
  points(x=x, y=c, pch=19, col="black")
  x=x+1
}
x=19.5
for(c in nb_2){
  points(x=x, y=c, pch=19, col="black")
  x=x+1
}
x=28.5
for(c in nb_3){
  points(x=x, y=c, pch=19, col="black")
  x=x+1
}

lines(x=c(1,9), y=c(70,70), col="red", lty=1, lwd=3)
lines(x=c(10,18), y=c(4,4), col="red", lty=1, lwd=3)
lines(x=c(19,27), y=c(2,2), col="red", lty=1, lwd=3)
lines(x=c(28,36), y=c(0,0), col="red", lty=1, lwd=3)
text(x=9.7, y=70, col="red", "70", cex=1.5)
text(x=18.7, y=4, col="red", "4", cex=1.5)
text(x=27.7, y=2, col="red", "2", cex=1.5)
text(x=36.7, y=0.65, col="red", "0", cex=1.5)

mtext("Number of genes", side=2, line=2.5, cex=1.8)
mtext("Number of LoF mutations per gene", side=1, line=3, cex=1.8)

#title(main="Distribution of the number of LoF mutations per gene",ylab="Number of genes", xlab="Number of LoF mutations per gene", cex.lab=1.8, cex.main=2)
#mtext("Number of genes", side=2, line=2, cex=2)




#legend(x=20, y=75, legend="Observed values", col="red", lty=1, bty="n", lwd=3, cex=1.5)

dev.off()


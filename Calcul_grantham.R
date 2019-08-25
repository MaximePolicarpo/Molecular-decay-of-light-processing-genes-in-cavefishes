library("boot")
library(plotrix)
library(data.table)
library(tidyverse)

read_plus <- function(flnm) {
  fread(file = flnm, fill = T, sep="\t") 
}

setwd("~/Desktop/Maximum_Likelihood_Mutpred_Pipeline/Results_Vision")


Astayanax_CF_mutpred <- list.files(path = ".", recursive = T, pattern = "Astyanax_CF_Vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))


Astayanax_SF_mutpred <- list.files(path = ".", recursive = T, pattern = "Astyanax_Surface_vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

Brotula_barbata_mutpred <- list.files(path = ".", recursive = T, pattern = "Brotula_barbata_vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

Carapus_acus_mutpred <- list.files(path = ".", recursive = T, pattern = "Carapus_acus_vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

Pygocentrus_mutpred <- list.files(path = ".", recursive = T, pattern = "Pygocentrus_nattereri_vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

Lamprogrammus_mutpred <- list.files(path = ".", recursive = T, pattern = "Lamprogrammus_exutus_vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

Ldentata_mutpred <- list.files(path = ".", recursive = T, pattern = "Lucifuga_dentata_vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

Lholguinensis_mutpred <- list.files(path = ".", recursive = T, pattern = "Lucifuga_holguinensis_vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))


anshuiensis_mutpred <- list.files(path = ".", recursive = T, pattern = "anshuiensis_vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

grahami_mutpred <- list.files(path = ".", recursive = T, pattern = "grahami_vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

rhinocerous_mutpred <- list.files(path = ".", recursive = T, pattern = "rhinocerous_vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

danio_rerio_mutpred <- list.files(path = ".", recursive = T, pattern = "Danio_rerio_vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

TEST_mutpred <- list.files(path = ".", recursive = T, pattern = "TEST.results", full.names = F) %>% # list files
  map_df(~read_plus(.))


setwd("~/Desktop/Maximum_Likelihood_Mutpred_Pipeline/SIMULATIONS/VISION")

Vision_Simulations_Danio <- list.files(path = ".", recursive = T, pattern = "Danio_rerio_Simulations_Vision.tsv", full.names = F) %>% # list files
  map_df(~read_plus(.))

Vision_Simulations_Pygocentrus <- list.files(path = ".", recursive = T, pattern = "Pygocentrus_Simulations_Vision.tsv", full.names = F) %>% # list files
  map_df(~read_plus(.))

Vision_Simulations_Brotula <- list.files(path = ".", recursive = T, pattern = "Brotula_barbata_Simulations_Vision.tsv", full.names = F) %>% # list files
  map_df(~read_plus(.))


colnames(Vision_Simulations_Danio) <- c("Simu", "gene", "Substitution", "MutPred Score")
colnames(Vision_Simulations_Pygocentrus) <- c("Simu", "gene", "Substitution", "MutPred Score")
colnames(Vision_Simulations_Brotula) <- c("Simu", "gene", "Substitution", "MutPred Score")



##### PLOT DE LA SIMULATION DANIO TS/TV : 1.56 ###

test_matrice <- matrix(nrow=300)
for(i in unique(Vision_Simulations_Danio$"Simu")){
  test_matrice <- cbind(test_matrice, sample(as.numeric(unlist(c(Vision_Simulations_Danio[Vision_Simulations_Danio$Simu==i][,4]))), 300))
}
test_matrice <- test_matrice[,-1]

simu_all <- c()
for(i in 1:100){
  simu_all <- c(simu_all, as.numeric(test_matrice[,i]))
}


test_matrice_sorted_col <- apply(test_matrice, 2, sort)
test_matrice_sorted <- t(apply(test_matrice_sorted_col, 1, sort))

Fmax <- ecdf(test_matrice_sorted[,100])
Fmin <- ecdf(test_matrice_sorted[,1])
F94min <- ecdf(test_matrice_sorted[,4])
F94max <- ecdf(test_matrice_sorted[,97])
F74min <- ecdf(test_matrice_sorted[,14])
F74max <- ecdf(test_matrice_sorted[,87])
F50min <- ecdf(test_matrice_sorted[,26])
F50max <- ecdf(test_matrice_sorted[,75])


plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,1), xlim=c(0,1), main="Empirical distribution of MutPred scores for vision genes mutations", cex.lab=1.5, cex.axis=1.7, cex.main=2.3)


for(i in 1:100){
  plot(ecdf(test_matrice[,i]), verticals = TRUE, do.points = FALSE, lwd=1, add=TRUE, col="gray93")
}


tmpx <- seq(min(test_matrice_sorted[,4],test_matrice_sorted[,97]), max(test_matrice_sorted[,4],test_matrice_sorted[,97]), len=10000) 
p0 <- F94min(tmpx)
p1 <- F94max(tmpx)
segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "gray88", "gray88"), alpha.f=0.1))

tmpx <- seq(min(test_matrice_sorted[,14],test_matrice_sorted[,87]), max(test_matrice_sorted[,14],test_matrice_sorted[,87]), len=10000) 
p0 <- F74min(tmpx)
p1 <- F74max(tmpx)
segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "gray78", "gray78"), alpha.f=0.1))


tmpx <- seq(min(test_matrice_sorted[,26],test_matrice_sorted[,75]), max(test_matrice_sorted[,26],test_matrice_sorted[,75]), len=10000) 
p0 <- F50min(tmpx)
p1 <- F50max(tmpx)
segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "gray69", "gray68"), alpha.f=0.1))



plot(ecdf(simu_all), verticals = TRUE, do.points = FALSE, lwd=2, add=TRUE)


## PLOT SIMULATION PYGO TS/TV : 2.3 ##

test_matrice <- matrix(nrow=58)
for(i in unique(Vision_Simulations_Pygocentrus$"Simu")){
  test_matrice <- cbind(test_matrice, sample(as.numeric(unlist(c(Vision_Simulations_Pygocentrus[Vision_Simulations_Pygocentrus$Simu==i][,4]))), 58))
}
test_matrice <- test_matrice[,-1]

simu_all <- c()
for(i in 1:100){
  simu_all <- c(simu_all, as.numeric(test_matrice[,i]))
}


test_matrice_sorted_col <- apply(test_matrice, 2, sort)
test_matrice_sorted <- t(apply(test_matrice_sorted_col, 1, sort))

Fmax <- ecdf(test_matrice_sorted[,100])
Fmin <- ecdf(test_matrice_sorted[,1])
F94min <- ecdf(test_matrice_sorted[,4])
F94max <- ecdf(test_matrice_sorted[,97])
F74min <- ecdf(test_matrice_sorted[,14])
F74max <- ecdf(test_matrice_sorted[,87])
F50min <- ecdf(test_matrice_sorted[,26])
F50max <- ecdf(test_matrice_sorted[,75])


for(i in 1:100){
  plot(ecdf(test_matrice[,i]), verticals = TRUE, do.points = FALSE, lwd=1, add=TRUE, col="lightblue")
}


tmpx <- seq(min(test_matrice_sorted[,4],test_matrice_sorted[,97]), max(test_matrice_sorted[,4],test_matrice_sorted[,97]), len=10000) 
p0 <- F94min(tmpx)
p1 <- F94max(tmpx)
segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "steelblue", "steelblue"), alpha.f=0.1))

tmpx <- seq(min(test_matrice_sorted[,14],test_matrice_sorted[,87]), max(test_matrice_sorted[,14],test_matrice_sorted[,87]), len=10000) 
p0 <- F74min(tmpx)
p1 <- F74max(tmpx)
segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "royalblue1", "royalblue1"), alpha.f=0.1))


tmpx <- seq(min(test_matrice_sorted[,26],test_matrice_sorted[,75]), max(test_matrice_sorted[,26],test_matrice_sorted[,75]), len=10000) 
p0 <- F50min(tmpx)
p1 <- F50max(tmpx)
segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "royalblue4", "royalblue4"), alpha.f=0.1))


plot(ecdf(simu_all), verticals = TRUE, do.points = FALSE, lwd=3, add=TRUE, col="goldenrod2")

## PLOT SIMULATION Brotula TS/TV : 4 ##

test_matrice <- matrix(nrow=255)
for(i in unique(Vision_Simulations_Brotula$"Simu")){
  test_matrice <- cbind(test_matrice, sample(as.numeric(unlist(c(Vision_Simulations_Brotula[Vision_Simulations_Brotula$Simu==i][,4]))), 255))
}
test_matrice <- test_matrice[,-1]

simu_all <- c()
for(i in 1:100){
  simu_all <- c(simu_all, as.numeric(test_matrice[,i]))
}


test_matrice_sorted_col <- apply(test_matrice, 2, sort)
test_matrice_sorted <- t(apply(test_matrice_sorted_col, 1, sort))

Fmax <- ecdf(test_matrice_sorted[,100])
Fmin <- ecdf(test_matrice_sorted[,1])
F94min <- ecdf(test_matrice_sorted[,4])
F94max <- ecdf(test_matrice_sorted[,97])
F74min <- ecdf(test_matrice_sorted[,14])
F74max <- ecdf(test_matrice_sorted[,87])
F50min <- ecdf(test_matrice_sorted[,26])
F50max <- ecdf(test_matrice_sorted[,75])


for(i in 1:100){
  plot(ecdf(test_matrice[,i]), verticals = TRUE, do.points = FALSE, lwd=1, add=TRUE, col="gold2")
}


tmpx <- seq(min(test_matrice_sorted[,4],test_matrice_sorted[,97]), max(test_matrice_sorted[,4],test_matrice_sorted[,97]), len=10000) 
p0 <- F94min(tmpx)
p1 <- F94max(tmpx)
segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "goldenrod3", "goldenrod3"), alpha.f=0.1))

tmpx <- seq(min(test_matrice_sorted[,14],test_matrice_sorted[,87]), max(test_matrice_sorted[,14],test_matrice_sorted[,87]), len=10000) 
p0 <- F74min(tmpx)
p1 <- F74max(tmpx)
segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "goldenrod4", "goldenrod4"), alpha.f=0.1))


tmpx <- seq(min(test_matrice_sorted[,26],test_matrice_sorted[,75]), max(test_matrice_sorted[,26],test_matrice_sorted[,75]), len=10000) 
p0 <- F50min(tmpx)
p1 <- F50max(tmpx)
segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "darkorange3", "darkorange3"), alpha.f=0.1))


plot(ecdf(simu_all), verticals = TRUE, do.points = FALSE, lwd=3, add=TRUE, col="deeppink3")


legend("topleft", title="Simulations", legend=c("D.rerio CDS (300) K=1.8", "P.nattereri CDS (54) K=2.3", "B.barbata CDS (255) K=4.5"), lty=1, col=c("black", "goldenrod2", "deeppink3"), bty="n")


#### PLOT DES VALEURS OBSERVEES ###

plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,1), xli=c(0,1),title.adj=c(0.4,1))


plot(ecdf(as.numeric(unlist(Astayanax_CF_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="firebrick1", lwd=3, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(Ldentata_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="brown4", lwd=3, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(Lholguinensis_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="orange", lwd=3, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(Astayanax_SF_mutpred[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="deeppink", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(Brotula_barbata_mutpred[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="chartreuse4", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(Lamprogrammus_mutpred[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="darkgreen", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(Carapus_acus_mutpred[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="aquamarine", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(Pygocentrus_mutpred[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="darkviolet", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(danio_rerio_mutpred[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="green", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(grahami_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="blue4", lwd=0.8, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(anshuiensis_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="darkred", lwd=0.8, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(rhinocerous_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="goldenrod1", lwd=0.8, cex=0.9, pch=2)


#ks.test(as.numeric(unlist(Astayanax_CF_mutpred[,3])), simu_all)
#ks.test(as.numeric(unlist(Astayanax_SF_mutpred[,3])), simu_all)
#ks.test(as.numeric(unlist(Ldentata_mutpred[,3])), simu_all)


legend("bottomright", legend=c("Simulations","Mean of 100 simulations", "Astayanax mexicanus SF (45)", "Pygocentrus nattereri (1658)","Danio rerio (1882)", "Brotula barbata (1406)", "Astyanax mexicanus CF (54)","Lucifuga gibarensis (120)", "Lucifuga dentata (280)"), col=c("gray93", "black", "deeppink", "darkviolet", "green", "chartreuse4", "firebrick1", "orange", "brown4"), lty=1, text.font=3, lwd=c(3, 3, 3, 3, 3, 3, 4, 4, 4), bty= 'n', pch=c(NA, NA, NA, NA, NA, NA, 17, 17, 17))
legend("bottomright", legend=c("Simulations","Mean of 100 simulations", "Sinocyclocheilus grahami (535)", "Sinocyclocheilus rhinocerous (612)","Sinocyclocheilus anshuiensis (553)"), col=c("gray93", "black", "blue4", "goldenrod1", "darkred"), lty=1, text.font=3, lwd=c(3, 3, 3, 3, 3), bty= 'n', pch=c(NA, NA, 17, 17, 17))


#legend("topleft", legend=c("Astyanax mexicanus CF","Lucifuga holguinensis", "Lucifuga dentata"), col=c("firebrick1", "orange", "brown4"), pch=2, text.font=3)
#legend("bottomright", legend=c("grahami", "rhinocerous", "anshuiensis"), col=c("blue4", "deepskyblue", "cyan"), bty="n", lty=1)

#legend("center", legend=c("Simulations", "94% interval", "74% interval", "50% interval", "Mean"), col=c("gray93", "gray88", "gray78", "gray69", "black"), pch=15)

dep1 <- rep("NUL", length(simu_all))

simu_all_matrix <- cbind.data.frame(dep1, dep1, simu_all)

n=0
matrix_pvalues <- matrix(nrow=8, ncol=8)
malist <- list(simu_all_matrix, Astayanax_SF_mutpred, Pygocentrus_mutpred, danio_rerio_mutpred, Brotula_barbata_mutpred, Astayanax_CF_mutpred, Lholguinensis_mutpred, Ldentata_mutpred)
tests_concat <- c()
for(i in malist){
  tests_concat <- c()
  for(j in malist){
    tests_concat <- c(tests_concat, ks.test(as.numeric(unlist(i[,3])), as.numeric(unlist(j[,3])))$p)
  }
  n=n+1
  matrix_pvalues[n,] <- tests_concat
}

#tests_concat <- round(tests_concat, digits=6)


matrix_pvalues <- round(matrix_pvalues, digits=6)

library(RColorBrewer)

nColors <- 1000
cols <-colorRampPalette(colors=c("forestgreen", "red"))(nColors)
zScale <- seq(min(matrix_pvalues), max(matrix_pvalues), length.out = nColors)

findNearestColour <- function(x) {
  colorIndex <- which(abs(zScale - x) == min(abs(zScale - x)))
  return(cols[colorIndex])
}

rows <- 8
collumns <- 8

plot(1, 1, type = "n", xlim = c(1, collumns), ylim = c(1, rows), 
     axes = F, xlab = "", ylab = "")

for(r in 1:rows){
  for(c in 1:collumns){
    text(c, r, matrix_pvalues[c, r], col = findNearestColour(matrix_pvalues[c, r])) 
  }
}


#cellcol<-color.scale(cbind(matrix_pvalues,c(-1,rep(1,7))),c(0.3,0,0.3,0.1,1,0.95,1))[,1:8]

cellcol<-matrix(rep("#000000",64),nrow=8)

#cellcol[matrix_pvalues<0.05]<-color.scale(matrix_pvalues[matrix_pvalues<0.05],c(1,0.8),c(0,0.8),0.1)
#cellcol[matrix_pvalues>0.05]<-color.scale(matrix_pvalues[matrix_pvalues>0.05],0,c(0.95,1),c(0.8,0))

cellcol[matrix_pvalues<0.05] <- "brown4"
cellcol[matrix_pvalues>0.05] <- "forestgreen"

color2D.matplot(matrix_pvalues, show.values=TRUE,cellcolors=cellcol, xlab="", ylab="",main="Kolmogorov-Smirnov Tests p-values between MutPred scores distribution", axes=FALSE, vcex=2.5, cex.main=2)



###################################################################################
###################### KERNEL DENSITIES #############################################
###################################################################################
###################################################################################




plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,3.2), xlim=c(0,1), cex.axis=2, cex.lab=1.5, main="Kernel densities of vision genes mutpred scores", cex.main=2.3)


matrice_de_y=matrix(nrow=512)
matrice_de_x=matrix(nrow=512)
mean_distrib_simu <- c()
for(i in unique(Vision_Simulations_Danio$"Simu")){
  matrice_de_x <- cbind(matrice_de_x, density(sample(as.numeric(unlist(c(Vision_Simulations_Danio[Vision_Simulations_Danio$Simu==i][,4]))), 300), from=0, to=1)$x)
  matrice_de_y <- cbind(matrice_de_y, density(sample(as.numeric(unlist(c(Vision_Simulations_Danio[Vision_Simulations_Danio$Simu==i][,4]))), 300), from=0, to=1)$y)
  mean_distrib_simu <- c(mean_distrib_simu, sample(as.numeric(unlist(c(Vision_Simulations_Danio[Vision_Simulations_Danio$Simu==i][,4]))), 300))
}


matrice_de_x <- matrice_de_x[,-1]
matrice_de_y <- matrice_de_y[,-1]
moyenne_kernel_x <- apply(matrice_de_x, 1, mean) 
moyenne_kernel_y <- apply(matrice_de_y, 1, mean) 
max_kernel_densities_x <- apply(matrice_de_x, 1, max) 
max_kernel_density_y <- apply(matrice_de_y, 1, max)
min_kernel_density_y <- apply(matrice_de_y, 1, min)

for(i in unique(Vision_Simulations_Danio$"Simu")){
  lines(density(sample(as.numeric(unlist(c(Vision_Simulations_Danio[Vision_Simulations_Danio$Simu==i][,4]))), 300), from=0, to=1), col="gray88")
}


#polygon(c(max_kernel_densities_x,rev(max_kernel_densities_x)),c(max_kernel_density_y,rev(min_kernel_density_y)), col="grey88", border=NA)



sorted_coordinates <- apply(matrice_de_y, 1, sort) 
min_values <- sorted_coordinates[11,]
max_values <- sorted_coordinates[90,]
polygon(c(max_kernel_densities_x,rev(max_kernel_densities_x)),c(max_values,rev(min_values)), col="grey70", border=NA)

lines(x=moyenne_kernel_x, y=moyenne_kernel_y, col="black", lwd=2)


lines(density(as.numeric(unlist(Astayanax_CF_mutpred[,3]))), col="firebrick1", lwd=3)
lines(density(as.numeric(unlist(Ldentata_mutpred[,3]))), col="brown4", lwd=3)
lines(density(as.numeric(unlist(Lholguinensis_mutpred[,3]))), col="orange", lwd=3)
lines(density(as.numeric(unlist(Astayanax_SF_mutpred[,3]))), col="deeppink", lwd=2)
lines(density(as.numeric(unlist(Brotula_barbata_mutpred[,3]))), col="chartreuse4", lwd=2)
lines(density(as.numeric(unlist(Lamprogrammus_mutpred[,3]))), col="darkgreen", lwd=2)
lines(density(as.numeric(unlist(Carapus_acus_mutpred[,3]))),  col="aquamarine", lwd=2)
lines(density(as.numeric(unlist(Pygocentrus_mutpred[,3]))), col="darkviolet", lwd=2)
lines(density(as.numeric(unlist(danio_rerio_mutpred[,3]))), col="green", lwd=2)
lines(density(as.numeric(unlist(grahami_mutpred[,3]))), col="blue4", lwd=2)
lines(density(as.numeric(unlist(anshuiensis_mutpred[,3]))), col="cyan", lwd=2)
lines(density(as.numeric(unlist(rhinocerous_mutpred[,3]))), col="deepskyblue", lwd=2)


legend("bottomright", legend=c("Simulations","Mean of 100 simulations", "Astayanax mexicanus SF (45)", "Pygocentrus nattereri (1658)","Danio rerio (1882)", "Brotula barbata (1406)", "Astyanax mexicanus CF (54)","Lucifuga gibarensis (120)", "Lucifuga dentata (280)"), col=c("gray93", "black", "deeppink", "darkviolet", "green", "chartreuse4", "firebrick1", "orange", "brown4"), lty=1, text.font=3, lwd=c(3, 3, 3, 3, 3, 3, 4, 4, 4), bty= 'n')

###################################################################################
###################### ADMITURE MODELS #############################################
###################################################################################
###################################################################################

#Ldentata_distrib <- as.numeric(unlist(Ldentata_mutpred[,3]))
#
#Positiv_mutations <- as.numeric(unlist(Brotula_barbata_mutpred[,3]))
#
#my_mutations <- sample(Positiv_mutations, 100, replace = FALSE, prob = NULL)
#
##plot(ecdf(my_mutations), verticals = TRUE, do.points = TRUE, add=TRUE, col="firebrick1", lwd=1.4, cex=0.9, pch=2)
#
#Neutral_mutations <- simu_all
#
#plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,1), xli=c(0,1))
#
#test_admixture <- c()
#for(i in 1:1000){
#  my_mutations_final <- my_mutations
#  #my_mutations <- sample(Positiv_mutations, 100, replace = FALSE, prob = NULL)
#  my_mutations_final <- c(my_mutations_final, sample(Neutral_mutations, i, replace = FALSE, prob = NULL))
#  #plot(ecdf(my_mutations), verticals = TRUE, do.points = FALSE, lwd=1, add=TRUE, col="black")
#  test_admixture <- c(test_admixture, ks.test(Ldentata_distrib, my_mutations_final)$p)
#  
#}
#
#x=1:1000
#
#plot(x, test_admixture, type="l")

#Ldentata_distrib <- as.numeric(unlist(Ldentata_mutpred[,3]))
#Positiv_mutations <- as.numeric(unlist(Brotula_barbata_mutpred[,3]))
#Neutral_mutations <- simu_all
#
#
#test_admixture <- matrix(nrow=100, ncol=14060)
#for(j in 1:100){
#  my_vector_test <- c()
#  for(i in 1:14060){
#    Positiv_mutations <- as.numeric(unlist(Brotula_barbata_mutpred[,3]))
#    my_mutations_final <- c(Positiv_mutations, sample(Neutral_mutations, i, replace = FALSE, prob = NULL))
#    my_vector_test <- c(my_vector_test, ks.test(Ldentata_distrib, my_mutations_final)$p)
#  }
#  test_admixture[j,] <- my_vector_test
#}
#
#
#
#
#moyenne_simulations <- apply(test_admixture, 2, mean, na.rm=TRUE)
#plot(moyenne_simulations, type="l")
#
#nb_2 <- which(moyenne_simulations==max(moyenne_simulations))
#nb_1 <- length(as.numeric(unlist(Brotula_barbata_mutpred[,3])))
#
#
#nb_de_mut_total <- nb_1 + nb_2
#
#ma_matrix_distribution <- matrix(nrow=100, ncol=nb_de_mut_total)
#my_vector_test <- c()
#for(j in 1:100){
#  Positiv_mutations <- as.numeric(unlist(Brotula_barbata_mutpred[,3]))
#  my_mutations_final <- c(Positiv_mutations, sample(Neutral_mutations, nb_2, replace = FALSE, prob = NULL))
#  ma_matrix_distribution[j,] <- my_mutations_final
#}
#
#
#plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,1), xli=c(0,1), main="Empirical distribution of MutPred scores")
#
#for(i in 1:100){
#  plot(ecdf(as.numeric(unlist(ma_matrix_distribution[i,]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="green", lwd=1.4, cex=0.4, pch=19)
#}
#plot(ecdf(as.numeric(unlist(Ldentata_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="brown4", lwd=1.4, cex=0.9, pch=2)
#
#legend("bottomright", legend=c("Simulations 1406 Brotula mutation + 2420 simulated neutral mutations", "Lucifuga dentata distribution"), col=c("green", "brown4"), lty=1)
#
####################################################################################
####################################################################################
####################################################################################
####################################################################################
#
#
#Lholguin_distrib <- as.numeric(unlist(Lholguinensis_mutpred[,3]))
#Positiv_mutations <- as.numeric(unlist(Brotula_barbata_mutpred[,3]))
#Neutral_mutations <- simu_all
#
#
#test_admixture2 <- matrix(nrow=100, ncol=14060)
#for(j in 1:100){
#  my_vector_test2 <- c()
#  for(i in 1:14060){
#    Positiv_mutations <- as.numeric(unlist(Brotula_barbata_mutpred[,3]))
#    my_mutations_final <- c(Positiv_mutations, sample(Neutral_mutations, i, replace = FALSE, prob = NULL))
#    my_vector_test2 <- c(my_vector_test2, ks.test(Lholguin_distrib, my_mutations_final)$p)
#  }
#  test_admixture2[j,] <- my_vector_test2
#}
#
#
#
#
#moyenne_simulations2 <- apply(test_admixture2, 2, mean, na.rm=TRUE)
#plot(moyenne_simulations2, type="l")
#
#nb_2 <- which(moyenne_simulations2==max(moyenne_simulations2))
#nb_1 <- length(as.numeric(unlist(Brotula_barbata_mutpred[,3])))
#
#
#nb_de_mut_total2 <- nb_1 + nb_2
#
#ma_matrix_distribution2 <- matrix(nrow=100, ncol=nb_de_mut_total2)
#my_vector_test <- c()
#for(j in 1:100){
#  Positiv_mutations <- as.numeric(unlist(Brotula_barbata_mutpred[,3]))
#  my_mutations_final <- c(Positiv_mutations, sample(Neutral_mutations, nb_2, replace = FALSE, prob = NULL))
#  ma_matrix_distribution2[j,] <- my_mutations_final
#}
#
#
#plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,1), xli=c(0,1), main="Empirical distribution of MutPred scores")
#
#for(i in 1:100){
#  plot(ecdf(as.numeric(unlist(ma_matrix_distribution2[i,]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="green", lwd=1.4, cex=0.4, pch=19)
#}
#plot(ecdf(as.numeric(unlist(Lholguinensis_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="brown4", lwd=1.4, cex=0.9, pch=2)
#
#
#
#
####################################################################################
####################################################################################
####################################################################################
####################################################################################
#
##Simulation avec un nb de mut constant = 54 
##Soit 54 mut neutres issus de la simulation soit x% de site neutres et 100-x% de sites 
## issus de A.mex SF
#
#Amex_CF_distrib <- as.numeric(unlist(Astayanax_CF_mutpred[,3]))
#Positiv_mutations <- as.numeric(unlist(Astayanax_SF_mutpred[,3]))
#Neutral_mutations <- simu_all
#
#test_admixture <- matrix(nrow=1000, ncol=46)
#for(j in 1:1000){
#  my_vector_test <- c()
#  for(i in 0:45){
#    Positiv_mutations <- as.numeric(unlist(Astayanax_SF_mutpred[,3]))
#    my_mutations_final <- c(sample(Positiv_mutations, 45-i, replace = FALSE, prob = NULL), sample(Neutral_mutations, i, replace = FALSE, prob = NULL))
#    my_vector_test <- c(my_vector_test, ks.test(Amex_CF_distrib, my_mutations_final)$p)
#  }
#  test_admixture[j,] <- my_vector_test
#}
#
#
#
#
#moyenne_simulations <- apply(test_admixture, 2, mean, na.rm=TRUE)
#plot(moyenne_simulations, type="l", xlab="Nombre de mutations neutres", ylab="Kolmogorov-test", main="Graphique A.mex (Ntotal=45 mutations)")
#
#nb_2 <- which(moyenne_simulations==max(moyenne_simulations))
#nb_1 <- length(Positiv_mutations)
#
#
#nb_de_mut_total <- 45
#
#ma_matrix_distribution <- matrix(nrow=1000, ncol=nb_de_mut_total)
#my_vector_test <- c()
#for(j in 1:1000){
#  Positiv_mutations <- as.numeric(unlist(Astayanax_SF_mutpred[,3]))
#  my_mutations_final <- c(sample(Positiv_mutations, 45-nb_2, replace = FALSE, prob = NULL), sample(Neutral_mutations, nb_2, replace = FALSE, prob = NULL))
#  ma_matrix_distribution[j,] <- my_mutations_final
#}
#
#
#plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,1), xli=c(0,1), main="Empirical distribution of MutPred scores")
#
#for(i in 1:1000){
#  plot(ecdf(as.numeric(unlist(ma_matrix_distribution[i,]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="green", lwd=1.4, cex=0.4, pch=19)
#}
#plot(ecdf(as.numeric(unlist(Astayanax_CF_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="brown4", lwd=1.4, cex=0.9, pch=2)
#
#
####################################################################################
####################################################################################
####################################################################################
####################################################################################
#
#Lholguin_distrib <- as.numeric(unlist(Lholguinensis_mutpred[,3]))
#Positiv_mutations <- as.numeric(unlist(Brotula_barbata_mutpred[,3]))
#Neutral_mutations <- simu_all
#
#
#test_admixture <- matrix(nrow=1000, ncol=1001)
#for(j in 1:1000){
#  my_vector_test <- c()
#  for(i in 0:1000){
#    Positiv_mutations <- as.numeric(unlist(Brotula_barbata_mutpred[,3]))
#    my_mutations_final <- c(sample(Positiv_mutations, 1000-i, replace = FALSE, prob = NULL), sample(Neutral_mutations, i, replace = FALSE, prob = NULL))
#    my_vector_test <- c(my_vector_test, ks.test(Lholguin_distrib, my_mutations_final)$p)
#  }
#  test_admixture[j,] <- my_vector_test
#}
#
#
#
#
#moyenne_simulations <- apply(test_admixture, 2, mean, na.rm=TRUE)
#lines(moyenne_simulations, col="blue")
#
#nb_2 <- which(moyenne_simulations==max(moyenne_simulations))
##nb_1 <- length(as.numeric(unlist(Brotula_barbata_mutpred[,3])))
#
#
#nb_de_mut_total <- 1000
#
#ma_matrix_distribution <- matrix(nrow=1000, ncol=nb_de_mut_total)
#my_vector_test <- c()
#for(j in 1:1000){
#  Positiv_mutations <- as.numeric(unlist(Brotula_barbata_mutpred[,3]))
#  my_mutations_final <- c(sample(Positiv_mutations, 1000-nb_2, replace = FALSE, prob = NULL), sample(Neutral_mutations, nb_2, replace = FALSE, prob = NULL))
#  ma_matrix_distribution[j,] <- my_mutations_final
#}
#
#
#plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,1), xli=c(0,1), main="Empirical distribution of MutPred scores")
#
#for(i in 1:1000){
#  plot(ecdf(as.numeric(unlist(ma_matrix_distribution[i,]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="green", lwd=1.4, cex=0.4, pch=19)
#}
#plot(ecdf(as.numeric(unlist(Lholguinensis_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="brown4", lwd=1.4, cex=0.9, pch=2)
#
####
#
#Ldentata_distrib <- as.numeric(unlist(Ldentata_mutpred[,3]))
#Positiv_mutations <- as.numeric(unlist(Brotula_barbata_mutpred[,3]))
#Neutral_mutations <- simu_all
#
#
#test_admixture <- matrix(nrow=1000, ncol=1001)
#for(j in 1:1000){
#  my_vector_test <- c()
#  for(i in 0:1000){
#    Positiv_mutations <- as.numeric(unlist(Brotula_barbata_mutpred[,3]))
#    my_mutations_final <- c(sample(Positiv_mutations, 1000-i, replace = FALSE, prob = NULL), sample(Neutral_mutations, i, replace = FALSE, prob = NULL))
#    my_vector_test <- c(my_vector_test, ks.test(Ldentata_distrib, my_mutations_final)$p)
#  }
#  test_admixture[j,] <- my_vector_test
#}
#
#
#moyenne_simulations <- apply(test_admixture, 2, mean, na.rm=TRUE)
#plot(moyenne_simulations, type="l", col="black")
#
#
#legend("topleft", legend=c("L.dent", "L.gibaren"), col=c("black", "blue"), lty=1)
#
#nb_2 <- which(moyenne_simulations==max(moyenne_simulations))
##nb_1 <- length(as.numeric(unlist(Brotula_barbata_mutpred[,3])))
#
#
#nb_de_mut_total <- 1000
#
#ma_matrix_distribution <- matrix(nrow=1000, ncol=nb_de_mut_total)
#my_vector_test <- c()
#for(j in 1:1000){
#  Positiv_mutations <- as.numeric(unlist(Brotula_barbata_mutpred[,3]))
#  my_mutations_final <- c(sample(Positiv_mutations, 1000-nb_2, replace = FALSE, prob = NULL), sample(Neutral_mutations, nb_2, replace = FALSE, prob = NULL))
#  ma_matrix_distribution[j,] <- my_mutations_final
#}
#
#
#plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,1), xli=c(0,1), main="Empirical distribution of MutPred scores")
#
#for(i in 1:1000){
#  plot(ecdf(as.numeric(unlist(ma_matrix_distribution[i,]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="forestgreen", lwd=1.4, cex=0.4, pch=19)
#}
#plot(ecdf(as.numeric(unlist(Ldentata_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="red", lwd=1.4, cex=0.9, pch=2)




Neutral_mutations <- simu_all
Positiv_mutations <- c(as.numeric(unlist(Brotula_barbata_mutpred[,3])), as.numeric(unlist(Pygocentrus_mutpred[,3])), as.numeric(unlist(danio_rerio_mutpred[,3])), as.numeric(unlist(Astayanax_SF_mutpred[,3])))
                       
Ldentata_distrib <- as.numeric(unlist(Ldentata_mutpred[,3]))
Amex_CF_distrib <- as.numeric(unlist(Astayanax_CF_mutpred[,3]))
Lholguin_distrib <- as.numeric(unlist(Lholguinensis_mutpred[,3]))


test_admixture_dent <- matrix(nrow=1000, ncol=1001)
test_admixture_holg <- matrix(nrow=1000, ncol=1001)
test_admixture_asty <- matrix(nrow=1000, ncol=1001)
for(j in 1:1000){
  my_vector_test_dent  <- c()
  my_vector_test_holg  <- c()
  my_vector_test_asty  <- c() 
  for(i in 0:1000){
    my_mutations_final <- c(sample(Positiv_mutations, 1000-i, replace = FALSE, prob = NULL), sample(Neutral_mutations, i, replace = FALSE, prob = NULL))
    my_vector_test_dent <- c(my_vector_test_dent, ks.test(Ldentata_distrib, my_mutations_final)$p)
    my_vector_test_holg <- c(my_vector_test_holg, ks.test(Lholguin_distrib, my_mutations_final)$p)
    my_vector_test_asty <- c(my_vector_test_asty, ks.test(Amex_CF_distrib, my_mutations_final)$p)
  }
  test_admixture_dent[j,] <- my_vector_test_dent
  test_admixture_holg[j,] <- my_vector_test_holg
  test_admixture_asty[j,] <- my_vector_test_asty
}


moyenne_simulations_dent <- apply(test_admixture_dent, 2, mean, na.rm=TRUE)
moyenne_simulations_holg <- apply(test_admixture_holg, 2, mean, na.rm=TRUE)
moyenne_simulations_asty <- apply(test_admixture_asty, 2, mean, na.rm=TRUE)

plot(moyenne_simulations_dent, type="l", col="black", xlab="Proportion of neutral mutations", ylab="Kolmogorov-test p-value", xaxt = "n", cex.lab=1.5, cex.axis=1.5, lwd=2, main="Kolmogorov-test p-value between observed cavefish distributions and an admixture of neutral and surface fish mutations")
lines(moyenne_simulations_holg, col="red", lwd=2)
lines(moyenne_simulations_asty, col="blue", lwd=2)
legend("topleft", legend=c("L.dentata", "L.gibarensis", "A.mexicanus CF"), col=c("black", "red", "blue"), lty=1, bty="n", lwd=2)

axis(1, at=c(0, 200, 400, 600, 800, 1000), labels=c("0", "0.2", "0.4", "0.6", "0.8", "1"), cex.axis=1.5)


nb_dent <- which(moyenne_simulations_dent==max(moyenne_simulations_dent))
nb_holg <- which(moyenne_simulations_holg==max(moyenne_simulations_holg))
nb_asty <- which(moyenne_simulations_asty==max(moyenne_simulations_asty))

ma_matrix_distribution_dent <- matrix(nrow=1000, ncol=1000)
ma_matrix_distribution_holg <- matrix(nrow=1000, ncol=1000)
ma_matrix_distribution_asty <- matrix(nrow=1000, ncol=1000)

for(j in 0:1000){
  Positiv_mutations <- as.numeric(unlist(Brotula_barbata_mutpred[,3]))
  my_mutations_final_dent <- c(sample(Positiv_mutations, 1000-nb_dent, replace = FALSE, prob = NULL), sample(Neutral_mutations, nb_dent, replace = FALSE, prob = NULL))
  my_mutations_final_holg <- c(sample(Positiv_mutations, 1000-nb_holg, replace = FALSE, prob = NULL), sample(Neutral_mutations, nb_holg, replace = FALSE, prob = NULL))
  my_mutations_final_asty <- c(sample(Positiv_mutations, 1000-nb_asty, replace = FALSE, prob = NULL), sample(Neutral_mutations, nb_asty, replace = FALSE, prob = NULL))
  
  ma_matrix_distribution_dent[j,] <- my_mutations_final_dent
  ma_matrix_distribution_holg[j,] <- my_mutations_final_holg
  ma_matrix_distribution_asty[j,] <- my_mutations_final_asty
}

plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,1), xlim=c(0,1), main="Empirical distribution of MutPred scores")

for(i in 1:100){
  plot(ecdf(as.numeric(unlist(ma_matrix_distribution_dent[i,]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="green", lwd=1.4, cex=0.4, pch=19)
  plot(ecdf(as.numeric(unlist(ma_matrix_distribution_holg[i,]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="darkseagreen1", lwd=1.4, cex=0.4, pch=19)
  plot(ecdf(as.numeric(unlist(ma_matrix_distribution_asty[i,]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="forestgreen", lwd=1.4, cex=0.4, pch=19)
}
plot(ecdf(as.numeric(unlist(Astayanax_CF_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="firebrick1", lwd=1.4, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(Ldentata_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="brown4", lwd=1.4, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(Lholguinensis_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="orange", lwd=1.4, cex=0.9, pch=2)




###################################################################################
###################################################################################
###################################################################################
###################################################################################


par(mfrow=c(1,3))

Grantham_matrix <- matrix(c(  0, 110, 145,  74,  58,  99, 124,  56, 142, 155, 144, 112,  89,  68,  46, 121,  65,  80, 135, 177,
                              110,   0, 102, 103,  71, 112,  96, 125,  97,  97,  77, 180,  29,  43,  86,  26,  96,  54,  91, 101,
                              145, 102,   0,  98,  92,  96,  32, 138,   5,  22,  36, 198,  99, 113, 153, 107, 172, 138,  15,  61,
                              74, 103,  98,   0,  38,  27,  68,  42,  95, 114, 110, 169,  77,  76,  91, 103, 108,  93,  87, 147,
                              58,  71,  92,  38,   0,  58,  69,  59,  89, 103,  92, 149,  47,  42,  65,  78,  85,  65,  81, 128,
                              99, 112,  96,  27,  58,   0,  64,  60,  94, 113, 112, 195,  86,  91, 111, 106, 126, 107,  84, 148,
                              124,  96,  32,  68,  69,  64,   0, 109,  29,  50,  55, 192,  84,  96, 133,  97, 152, 121,  21,  88,
                              56, 125, 138,  42,  59,  60, 109,   0, 135, 153, 147, 159,  98,  87,  80, 127,  94,  98, 127, 184,
                              142,  97,   5,  95,  89,  94,  29, 135,   0,  21,  33, 198,  94, 109, 149, 102, 168, 134,  10,  61,
                              155,  97,  22, 114, 103, 113,  50, 153,  21,   0,  22, 205, 100, 116, 158, 102, 177, 140,  28,  40,
                              144,  77,  36, 110,  92, 112,  55, 147,  33,  22,   0, 194,  83,  99, 143,  85, 160, 122,  36,  37,
                              112, 180, 198, 169, 149, 195, 192, 159, 198, 205, 194,   0, 174, 154, 139, 202, 154, 170, 196, 215,
                              89,  29,  99,  77,  47,  86,  84,  98,  94, 100,  83, 174,   0,  24,  68,  32,  81,  40,  87, 115,
                              68,  43, 113,  76,  42,  91,  96,  87, 109, 116,  99, 154,  24,   0,  46,  53,  61,  29, 101, 130,
                              46,  86, 153,  91,  65, 111, 133,  80, 149, 158, 143, 139,  68,  46,   0,  94,  23,  42, 142, 174,
                              121,  26, 107, 103,  78, 106,  97, 127, 102, 102,  85, 202,  32,  53,  94,   0, 101,  56,  95, 110,
                              65,  96, 172, 108,  85, 126, 152,  94, 168, 177, 160, 154,  81,  61,  23, 101,   0,  45, 160, 181,
                              80,  54, 138,  93,  65, 107, 121,  98, 134, 140, 122, 170,  40,  29,  42,  56,  45,   0, 126, 152,
                              135,  91,  15,  87,  81,  84,  21, 127,  10,  28,  36, 196,  87, 101, 142,  95, 160, 126,   0,  67,
                              177, 101,  61, 147, 128, 148,  88, 184,  61,  40,  37, 215, 115, 130, 174, 110, 181, 152,  67,   0), 
                          nrow=20, byrow=TRUE)
rownames(Grantham_matrix) <- c("S", "R", "L", "P", "T", "A", "V", "G", "I", "F", "Y", "C", "H", "Q", "N", "K", "D", "E", "M", "W")
colnames(Grantham_matrix) <- rownames(Grantham_matrix)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "Grantham Score", ylim=c(0,1), xli=c(0,250))

simu_all_grantham <- c()
for(i in unique(Vision_Simulations_Danio$"Simu")){
  aa1_list <- c()
  aa2_list <- c()
  for(j in unlist(Vision_Simulations_Danio[Vision_Simulations_Danio$Simu==i][,3])){
    aa1_list <- c(aa1_list, substr(j, 1, 1))
    aa2_list <- c(aa2_list, substrRight(j, 1))
  }
  grantham_list <- c()
  for(x in 1:length(aa1_list)){
    grantham_list <- c(grantham_list, Grantham_matrix[aa1_list[x], aa2_list[x]])
  }
  #plot(ecdf(grantham_list), verticals = TRUE, do.points = TRUE, add=TRUE, col="grey", lwd=1.4, cex=0.4, pch=19)
  #lines(density(grantham_list, from=-30, to=250), col="gray", lwd=2, add=TRUE)
  simu_all_grantham <- c(simu_all_grantham, grantham_list)
}

plot(density(grantham_list, from=-30, to=250), col="gray", lwd=2, ylim=c(0,0.02))


par(mfrow=c(1,1))
plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,2.5), xli=c(0,1))

lines(density(simu_all, from=0, to=1), col="black", lwd=2)




malist <- list(simu_all_grantham, grantham_list)

multhist(malist, freq=FALSE, col=c("gray", "forestgreen"), xlab="Grantham scores", ylab="Frequency", ylim=c(0,0.020))






###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(Astayanax_CF_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

grantham_list <- c()
for(i in 1:length(aa1_list)){
  grantham_list <- c(grantham_list, Grantham_matrix[aa1_list[i], aa2_list[i]])
}

plot(ecdf(grantham_list), verticals = TRUE, do.points = TRUE, add=TRUE, col="firebrick1", lwd=1.4, cex=0.4, pch=19)





###################################################################################

aa1_list <- c()
aa2_list <- c()

for(i in unlist(Astayanax_SF_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

grantham_list <- c()
for(i in 1:length(aa1_list)){
  grantham_list <- c(grantham_list, Grantham_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(grantham_list), verticals = TRUE, do.points = TRUE, add=TRUE, col="deeppink", lwd=1.4, cex=0.4, pch=19)


###################################################################################

aa1_list <- c()
aa2_list <- c()

for(i in unlist(Ldentata_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

grantham_list <- c()
for(i in 1:length(aa1_list)){
  grantham_list <- c(grantham_list, Grantham_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(grantham_list), verticals = TRUE, do.points = TRUE, add=TRUE, col="brown4", lwd=1.4, cex=0.4, pch=19)

###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(Lholguinensis_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

grantham_list <- c()
for(i in 1:length(aa1_list)){
  grantham_list <- c(grantham_list, Grantham_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(grantham_list), verticals = TRUE, do.points = TRUE, add=TRUE, col="orange", lwd=1.4, cex=0.4, pch=19)



###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(Brotula_barbata_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

grantham_list <- c()
for(i in 1:length(aa1_list)){
  grantham_list <- c(grantham_list, Grantham_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(grantham_list), verticals = TRUE, do.points = TRUE, add=TRUE, col="green", lwd=1.4, cex=0.4, pch=19)


###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(Lamprogrammus_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

grantham_list <- c()
for(i in 1:length(aa1_list)){
  grantham_list <- c(grantham_list, Grantham_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(grantham_list), verticals = TRUE, do.points = TRUE, add=TRUE, col="chartreuse4", lwd=1.4, cex=0.4, pch=19)


###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(Carapus_acus_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

grantham_list <- c()
for(i in 1:length(aa1_list)){
  grantham_list <- c(grantham_list, Grantham_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(grantham_list), verticals = TRUE, do.points = TRUE, add=TRUE, col="darkgreen", lwd=1.4, cex=0.4, pch=19)

###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(Pygocentrus_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

grantham_list <- c()
for(i in 1:length(aa1_list)){
  grantham_list <- c(grantham_list, Grantham_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(grantham_list), verticals = TRUE, do.points = TRUE, add=TRUE, col="aquamarine", lwd=1.4, cex=0.4, pch=19)



###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(danio_rerio_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

grantham_list <- c()
for(i in 1:length(aa1_list)){
  grantham_list <- c(grantham_list, Grantham_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(grantham_list), verticals = TRUE, do.points = TRUE, add=TRUE, col="forestgreen", lwd=1.4, cex=0.4, pch=19)


#lines(density(grantham_list, from=-30, to=250), col="forestgreen", lwd=2, add=TRUE)

###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(grahami_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

grantham_list <- c()
for(i in 1:length(aa1_list)){
  grantham_list <- c(grantham_list, Grantham_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(grantham_list), verticals = TRUE, do.points = TRUE, add=TRUE, col="green", lwd=1.4, cex=0.4, pch=19)


###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(anshuiensis_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

grantham_list <- c()
for(i in 1:length(aa1_list)){
  grantham_list <- c(grantham_list, Grantham_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(grantham_list), verticals = TRUE, do.points = TRUE, add=TRUE, col="blue4", lwd=1.4, cex=0.4, pch=19)


###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(rhinocerous_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

grantham_list <- c()
for(i in 1:length(aa1_list)){
  grantham_list <- c(grantham_list, Grantham_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(grantham_list), verticals = TRUE, do.points = TRUE, add=TRUE, col="deepskyblue", lwd=1.4, cex=0.4, pch=19)




legend("bottomright", col=c("grey", "firebrick1", "deeppink", "brown4","orange","green", "aquamarine","darkviolet"), legend=c("Simulations", "A.mex CF", "A.mex SF", "L.dent", "L.gib", "B. barbata", "P. nat", "D.rerio"), lty=1)
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################

setwd("~/Desktop/Maximum_Likelihood_Mutpred_Pipeline/Results_Circadian")


Astayanax_CF_mutpred_circ <- list.files(path = ".", recursive = T, pattern = "Astyanax_mexicanus_Cave_circadian.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))


Astayanax_SF_mutpred_circ <- list.files(path = ".", recursive = T, pattern = "Astyanax_mexicanus_Surface_circadian.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

Brotula_barbata_mutpred_circ <- list.files(path = ".", recursive = T, pattern = "Brotula_barbata_circadian.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

Carapus_acus_mutpred_circ <- list.files(path = ".", recursive = T, pattern = "Carapus_acus_circadian.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

Pygocentrus_mutpred_circ <- list.files(path = ".", recursive = T, pattern = "Pygocentrus_nattereri_circadian.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

Lamprogrammus_mutpred_circ <- list.files(path = ".", recursive = T, pattern = "Lamprogrammus_exutus_circadian.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

Ldentata_mutpred_circ <- list.files(path = ".", recursive = T, pattern = "Lucifuga_dentata_circadian.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

Lholguinensis_mutpred_circ <- list.files(path = ".", recursive = T, pattern = "Lucifuga_holguinensis_circadian.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))


anshuiensis_mutpred_circ <- list.files(path = ".", recursive = T, pattern = "Sinocyclocheilus_anshuiensis_circadian.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

grahami_mutpred_circ <- list.files(path = ".", recursive = T, pattern = "Sinocyclocheilus_grahami_circadian.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

rhinocerous_mutpred_circ <- list.files(path = ".", recursive = T, pattern = "Sinocyclocheilus_rhinocerous_circadian.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

danio_rerio_mutpred_circ <- list.files(path = ".", recursive = T, pattern = "Danio_rerio_circadian.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

setwd("~/Desktop/Maximum_Likelihood_Mutpred_Pipeline/SIMULATIONS/CIRCADIAN")

Circadian_Simulations_Danio <- list.files(path = ".", recursive = T, pattern = "Danio_rerio_Simulations_Circadian.tsv", full.names = F) %>% # list files
  map_df(~read_plus(.))

colnames(Circadian_Simulations_Danio) <- c("Simu", "gene", "Substitution", "MutPred Score")



####

test_matrice <- matrix(nrow=200)
for(i in unique(Circadian_Simulations_Danio$"Simu")){
  test_matrice <- cbind(test_matrice, sample(as.numeric(unlist(c(Circadian_Simulations_Danio[Circadian_Simulations_Danio$Simu==i][,4]))), 200))
}
test_matrice <- test_matrice[,-1]

simu_all <- c()
for(i in 1:100){
  simu_all <- c(simu_all, as.numeric(test_matrice[,i]))
}


test_matrice_sorted_col <- apply(test_matrice, 2, sort)
test_matrice_sorted <- t(apply(test_matrice_sorted_col, 1, sort))

Fmax <- ecdf(test_matrice_sorted[,100])
Fmin <- ecdf(test_matrice_sorted[,1])
F94min <- ecdf(test_matrice_sorted[,4])
F94max <- ecdf(test_matrice_sorted[,97])
F74min <- ecdf(test_matrice_sorted[,14])
F74max <- ecdf(test_matrice_sorted[,87])
F50min <- ecdf(test_matrice_sorted[,26])
F50max <- ecdf(test_matrice_sorted[,75])


plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,1), xli=c(0,1), main="Empirical distribution of MutPred scores for circadian genes mutations", cex.lab=1.5, cex.axis=1.7, cex.main=2.3)


for(i in 1:100){
  plot(ecdf(test_matrice[,i]), verticals = TRUE, do.points = FALSE, lwd=1, add=TRUE, col="gray93")
}


tmpx <- seq(min(test_matrice_sorted[,4],test_matrice_sorted[,97]), max(test_matrice_sorted[,4],test_matrice_sorted[,97]), len=10000) 
p0 <- F94min(tmpx)
p1 <- F94max(tmpx)
segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "gray88", "gray88"), alpha.f=0.1))

tmpx <- seq(min(test_matrice_sorted[,14],test_matrice_sorted[,87]), max(test_matrice_sorted[,14],test_matrice_sorted[,87]), len=10000) 
p0 <- F74min(tmpx)
p1 <- F74max(tmpx)
segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "gray78", "gray78"), alpha.f=0.1))


tmpx <- seq(min(test_matrice_sorted[,26],test_matrice_sorted[,75]), max(test_matrice_sorted[,26],test_matrice_sorted[,75]), len=10000) 
p0 <- F50min(tmpx)
p1 <- F50max(tmpx)
segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "gray69", "gray68"), alpha.f=0.1))


plot(ecdf(simu_all), verticals = TRUE, do.points = FALSE, lwd=2, add=TRUE)


###






plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,1), xli=c(0,1))


plot(ecdf(as.numeric(unlist(Astayanax_CF_mutpred_circ[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="firebrick1", lwd=3, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(Ldentata_mutpred_circ[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="brown4", lwd=3, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(Lholguinensis_mutpred_circ[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="orange", lwd=3, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(Astayanax_SF_mutpred_circ[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="deeppink", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(Brotula_barbata_mutpred_circ[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="chartreuse4", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(Lamprogrammus_mutpred_circ[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="darkgreen", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(Carapus_acus_mutpred_circ[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="aquamarine", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(Pygocentrus_mutpred_circ[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="darkviolet", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(danio_rerio_mutpred_circ[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="green", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(grahami_mutpred_circ[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="blue4", lwd=1.4, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(anshuiensis_mutpred_circ[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="darkred", lwd=1.4, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(rhinocerous_mutpred_circ[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="goldenrod1", lwd=1.4, cex=0.9, pch=2)

legend("bottomright", legend=c("Simulations","Mean of 100 simulations", "Astayanax mexicanus SF (38)", "Pygocentrus nattereri (888)","Danio rerio (1305)", "Brotula barbata (1234)", "Astyanax mexicanus CF (36)","Lucifuga gibarensis (56)", "Lucifuga dentata (152)"), col=c("gray93", "black", "deeppink", "darkviolet", "green", "chartreuse4", "firebrick1", "orange", "brown4"), lty=1, text.font=3, lwd=c(3, 3, 3, 3, 3, 3, 4, 4, 4), bty= 'n', pch=c(NA, NA, NA, NA, NA, NA, 17, 17, 17))
legend("bottomright", legend=c("Simulations","Mean of 100 simulations", "Sinocyclocheilus grahami (594)", "Sinocyclocheilus rhinocerous (594)","Sinocyclocheilus anshuiensis (576)"), col=c("gray93", "black", "blue4", "goldenrod1", "darkred"), lty=1, text.font=3, lwd=c(3, 3, 3, 3, 3), bty= 'n', pch=c(NA, NA, 17, 17, 17))



dep1 <- rep("NUL", length(simu_all))

simu_all_matrix <- cbind.data.frame(dep1, dep1, simu_all)


n=0
matrix_pvalues <- matrix(nrow=8, ncol=8)
malist <- list(simu_all_matrix, Astayanax_SF_mutpred_circ, Pygocentrus_mutpred_circ, danio_rerio_mutpred_circ, Brotula_barbata_mutpred_circ, Astayanax_CF_mutpred_circ, Lholguinensis_mutpred_circ, Ldentata_mutpred_circ)
tests_concat <- c()
for(i in malist){
  tests_concat <- c()
  for(j in malist){
    tests_concat <- c(tests_concat, ks.test(as.numeric(unlist(i[,3])), as.numeric(unlist(j[,3])))$p)
  }
  n=n+1
  matrix_pvalues[n,] <- tests_concat
}

cellcol[matrix_pvalues<0.05] <- "brown4"
cellcol[matrix_pvalues>0.05] <- "forestgreen"

color2D.matplot(matrix_pvalues, show.values=TRUE,cellcolors=cellcol, xlab="", ylab="",main="Kolmogorov-Smirnov Tests p-values between MutPred scores distribution", axes=FALSE, vcex=2.5, cex.main=2)




#legend("bottomright", legend=c("grahami", "rhinocerous", "anshuiensis"), col=c("blue4", "deepskyblue", "cyan"), bty="n", lty=1)

plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,1), xli=c(0,1))


Astayanax_CF_mutpred_circ_cry_per <- c(as.numeric(unlist(Astayanax_CF_mutpred_circ[Astayanax_CF_mutpred_circ$ID %like% "cry", ][,3])), as.numeric(unlist(Astayanax_CF_mutpred_circ[Astayanax_CF_mutpred_circ$ID %like% "per", ][,3])))
Astayanax_SF_mutpred_circ_cry_per <- c(as.numeric(unlist(Astayanax_SF_mutpred_circ[Astayanax_SF_mutpred_circ$ID %like% "cry", ][,3])), as.numeric(unlist(Astayanax_SF_mutpred_circ[Astayanax_SF_mutpred_circ$ID %like% "per", ][,3])))
Danio_circ_cry_per <- c(as.numeric(unlist(danio_rerio_mutpred_circ[danio_rerio_mutpred_circ$ID %like% "cry", ][,3])), as.numeric(unlist(danio_rerio_mutpred_circ[danio_rerio_mutpred_circ$ID %like% "per", ][,3])))
Brotula_circ_cry_per <- c(as.numeric(unlist(Brotula_barbata_mutpred_circ[Brotula_barbata_mutpred_circ$ID %like% "cry", ][,3])), as.numeric(unlist(Brotula_barbata_mutpred_circ[Brotula_barbata_mutpred_circ$ID %like% "per", ][,3])))
Pygocentrus_circ_cry_per <- c(as.numeric(unlist(Pygocentrus_mutpred_circ[Pygocentrus_mutpred_circ$ID %like% "cry", ][,3])), as.numeric(unlist(Pygocentrus_mutpred_circ[Pygocentrus_mutpred_circ$ID %like% "per", ][,3])))
dentata_circ_cry_per <- c(as.numeric(unlist(Ldentata_mutpred_circ[Ldentata_mutpred_circ$ID %like% "cry", ][,3])), as.numeric(unlist(Ldentata_mutpred_circ[Ldentata_mutpred_circ$ID %like% "per", ][,3])))
holguin_circ_cry_per <- c(as.numeric(unlist(Lholguinensis_mutpred_circ[Lholguinensis_mutpred_circ$ID %like% "cry", ][,3])), as.numeric(unlist(Lholguinensis_mutpred_circ[Lholguinensis_mutpred_circ$ID %like% "per", ][,3])))
simu_all_cry_per <- c(as.numeric(unlist(Circadian_Simulations_Danio[Circadian_Simulations_Danio$gene %like% "cry", ][,4])), as.numeric(unlist(Circadian_Simulations_Danio[Circadian_Simulations_Danio$gene %like% "per", ][,4])))


plot(ecdf(simu_all_cry_per), verticals = TRUE, do.points = TRUE, add=TRUE, col="black", lwd=1.4, cex=0.9, pch=2)
plot(ecdf(Astayanax_CF_mutpred_circ_cry_per), verticals = TRUE, do.points = TRUE, add=TRUE, col="firebrick1", lwd=1.4, cex=0.9, pch=2)
plot(ecdf(Astayanax_SF_mutpred_circ_cry_per), verticals = TRUE, do.points = TRUE, add=TRUE, col="deeppink", lwd=1.4, cex=0.9, pch=2)
plot(ecdf(Danio_circ_cry_per), verticals = TRUE, do.points = TRUE, add=TRUE, col="green", lwd=1.4, cex=0.9, pch=2)
plot(ecdf(Brotula_circ_cry_per), verticals = TRUE, do.points = TRUE, add=TRUE, col="chartreuse4", lwd=1.4, cex=0.9, pch=2)
plot(ecdf(Pygocentrus_circ_cry_per), verticals = TRUE, do.points = TRUE, add=TRUE, col="darkviolet", lwd=1.4, cex=0.9, pch=2)
plot(ecdf(dentata_circ_cry_per), verticals = TRUE, do.points = TRUE, add=TRUE, col="brown4", lwd=1.4, cex=0.9, pch=2)
plot(ecdf(holguin_circ_cry_per), verticals = TRUE, do.points = TRUE, add=TRUE, col="orange", lwd=1.4, cex=0.9, pch=2)


plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,3.2), xli=c(0,1))
lines(density(simu_all_cry_per), col="black", lwd=2)
lines(density(Astayanax_CF_mutpred_circ_cry_per), col="firebrick1", lwd=2)
lines(density(Astayanax_SF_mutpred_circ_cry_per), col="deeppink", lwd=2)



ks.test(Astayanax_CF_mutpred_circ_cry_per, simu_all_cry_per)
ks.test(as.numeric(unlist(Astayanax_CF_mutpred_circ[,3])), simu_all)

#### KERNEL DENSITIES ##

plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,8), xlim=c(0,1), cex.axis=2, cex.lab=1.5, main="Kernel densities of circadian clock genes mutpred scores", cex.main=2.3)


matrice_de_y=matrix(nrow=512)
matrice_de_x=matrix(nrow=512)
mean_distrib_simu <- c()
for(i in unique(Circadian_Simulations_Danio$"Simu")){
  matrice_de_x <- cbind(matrice_de_x, density(sample(as.numeric(unlist(c(Circadian_Simulations_Danio[Circadian_Simulations_Danio$Simu==i][,4]))), 200), from=0, to=1)$x)
  matrice_de_y <- cbind(matrice_de_y, density(sample(as.numeric(unlist(c(Circadian_Simulations_Danio[Circadian_Simulations_Danio$Simu==i][,4]))), 200), from=0, to=1)$y)
  mean_distrib_simu <- c(mean_distrib_simu, sample(as.numeric(unlist(c(Circadian_Simulations_Danio[Circadian_Simulations_Danio$Simu==i][,4]))), 200))
}


matrice_de_x <- matrice_de_x[,-1]
matrice_de_y <- matrice_de_y[,-1]
moyenne_kernel_x <- apply(matrice_de_x, 1, mean) 
moyenne_kernel_y <- apply(matrice_de_y, 1, mean) 
max_kernel_densities_x <- apply(matrice_de_x, 1, max) 
max_kernel_density_y <- apply(matrice_de_y, 1, max)
min_kernel_density_y <- apply(matrice_de_y, 1, min)

for(i in unique(Circadian_Simulations_Danio$"Simu")){
  lines(density(sample(as.numeric(unlist(c(Circadian_Simulations_Danio[Circadian_Simulations_Danio$Simu==i][,4]))), 200), from=0, to=1), col="gray88")
}


#polygon(c(max_kernel_densities_x,rev(max_kernel_densities_x)),c(max_kernel_density_y,rev(min_kernel_density_y)), col="grey88", border=NA)



sorted_coordinates <- apply(matrice_de_y, 1, sort) 
min_values <- sorted_coordinates[11,]
max_values <- sorted_coordinates[90,]
polygon(c(max_kernel_densities_x,rev(max_kernel_densities_x)),c(max_values,rev(min_values)), col="grey70", border=NA)

lines(x=moyenne_kernel_x, y=moyenne_kernel_y, col="black", lwd=2)


lines(density(as.numeric(unlist(Astayanax_CF_mutpred_circ[,3]))), col="firebrick1", lwd=3)
lines(density(as.numeric(unlist(Ldentata_mutpred_circ[,3]))), col="brown4", lwd=3)
lines(density(as.numeric(unlist(Lholguinensis_mutpred_circ[,3]))), col="orange", lwd=3)
lines(density(as.numeric(unlist(Astayanax_SF_mutpred_circ[,3]))), col="deeppink", lwd=2)
lines(density(as.numeric(unlist(Brotula_barbata_mutpred_circ[,3]))), col="chartreuse4", lwd=2)
lines(density(as.numeric(unlist(Lamprogrammus_mutpred_circ[,3]))), col="darkgreen", lwd=2)
lines(density(as.numeric(unlist(Carapus_acus_mutpred_circ[,3]))),  col="aquamarine", lwd=2)
lines(density(as.numeric(unlist(Pygocentrus_mutpred_circ[,3]))), col="darkviolet", lwd=2)
lines(density(as.numeric(unlist(danio_rerio_mutpred_circ[,3]))), col="green", lwd=2)
lines(density(as.numeric(unlist(grahami_mutpred_circ[,3]))), col="blue4", lwd=2)
lines(density(as.numeric(unlist(anshuiensis_mutpred_circ[,3]))), col="cyan", lwd=2)
lines(density(as.numeric(unlist(rhinocerous_mutpred_circ[,3]))), col="deepskyblue", lwd=2)


legend("bottomright", legend=c("Simulations","Mean of 100 simulations", "Astayanax mexicanus SF (38)", "Pygocentrus nattereri (888)","Danio rerio (1305)", "Brotula barbata (1234)", "Astyanax mexicanus CF (36)","Lucifuga gibarensis (56)", "Lucifuga dentata (152)"), col=c("gray93", "black", "deeppink", "darkviolet", "green", "chartreuse4", "firebrick1", "orange", "brown4"), lty=1, text.font=3, lwd=c(3, 3, 3, 3, 3, 3, 4, 4, 4), bty= 'n')








###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################


setwd("~/Desktop/Maximum_Likelihood_Mutpred_Pipeline/Results_Pigmentation")

Astayanax_CF_mutpred_pigm <- list.files(path = ".", recursive = T, pattern = "Astyanax_mexicanus_Cave_pigment.results.parsed_results", full.names = F) %>% # list files
  map_df(~read_plus(.))


Astayanax_SF_mutpred_pigm <- list.files(path = ".", recursive = T, pattern = "Astyanax_mexicanus_Surface_pigment.results.parsed_results", full.names = F) %>% # list files
  map_df(~read_plus(.))

Brotula_barbata_mutpred_pigm <- list.files(path = ".", recursive = T, pattern = "Brotula_barbata_pigment.results.parsed_results", full.names = F) %>% # list files
  map_df(~read_plus(.))

Carapus_acus_mutpred_pigm <- list.files(path = ".", recursive = T, pattern = "Carapus_acus_pigment.results.parsed_results", full.names = F) %>% # list files
  map_df(~read_plus(.))

Pygocentrus_mutpred_pigm <- list.files(path = ".", recursive = T, pattern = "Pygocentrus_pigment.results.parsed_results", full.names = F) %>% # list files
  map_df(~read_plus(.))

Lamprogrammus_mutpred_pigm <- list.files(path = ".", recursive = T, pattern = "Lamprogrammus_pigment.results.parsed_results", full.names = F) %>% # list files
  map_df(~read_plus(.))

Ldentata_mutpred_pigm <- list.files(path = ".", recursive = T, pattern = "Lucifuga_dentata_pigment.results.parsed_results", full.names = F) %>% # list files
  map_df(~read_plus(.))

Lholguinensis_mutpred_pigm <- list.files(path = ".", recursive = T, pattern = "Lucifuga_holguinensis_pigment.results.parsed_results", full.names = F) %>% # list files
  map_df(~read_plus(.))


anshuiensis_mutpred_pigm <- list.files(path = ".", recursive = T, pattern = "Sinocyclocheilus_anshuiensis_pigment.results.parsed_results", full.names = F) %>% # list files
  map_df(~read_plus(.))

grahami_mutpred_pigm <- list.files(path = ".", recursive = T, pattern = "Sinocyclocheilus_grahami_pigment.results.parsed_results", full.names = F) %>% # list files
  map_df(~read_plus(.))

rhinocerous_mutpred_pigm <- list.files(path = ".", recursive = T, pattern = "Sinocyclocheilus_rhinocerous_pigment.results.parsed_results", full.names = F) %>% # list files
  map_df(~read_plus(.))

danio_rerio_mutpred_pigm <- list.files(path = ".", recursive = T, pattern = "Danio_rerio_pigment.results.parsed_results", full.names = F) %>% # list files
  map_df(~read_plus(.))



setwd("~/Desktop/Maximum_Likelihood_Mutpred_Pipeline/SIMULATIONS/PIGMENTATION")

Pigmentation_Simulations_Danio <- list.files(path = ".", recursive = T, pattern = "Danio_rerio_Simulations_Pigmentation.tsv", full.names = F) %>% # list files
  map_df(~read_plus(.))

colnames(Pigmentation_Simulations_Danio) <- c("Simu", "gene", "Substitution", "MutPred Score")

###############

test_matrice <- matrix(nrow=1000)
for(i in unique(Pigmentation_Simulations_Danio$"Simu")){
  test_matrice <- cbind(test_matrice, sample(as.numeric(unlist(c(Pigmentation_Simulations_Danio[Pigmentation_Simulations_Danio$Simu==i][,4]))), 1000))
}
test_matrice <- test_matrice[,-1]

simu_all <- c()
for(i in 1:100){
  simu_all <- c(simu_all, as.numeric(test_matrice[,i]))
}


test_matrice_sorted_col <- apply(test_matrice, 2, sort)
test_matrice_sorted <- t(apply(test_matrice_sorted_col, 1, sort))

Fmax <- ecdf(test_matrice_sorted[,100])
Fmin <- ecdf(test_matrice_sorted[,1])
F94min <- ecdf(test_matrice_sorted[,4])
F94max <- ecdf(test_matrice_sorted[,97])
F74min <- ecdf(test_matrice_sorted[,14])
F74max <- ecdf(test_matrice_sorted[,87])
F50min <- ecdf(test_matrice_sorted[,26])
F50max <- ecdf(test_matrice_sorted[,75])


plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,1), xli=c(0,1), main="Empirical distribution of MutPred scores for pigmentation genes mutations", cex.lab=1.5, cex.axis=1.7, cex.main=2.1)


for(i in 1:100){
  plot(ecdf(test_matrice[,i]), verticals = TRUE, do.points = FALSE, lwd=1, add=TRUE, col="gray93")
}


tmpx <- seq(min(test_matrice_sorted[,4],test_matrice_sorted[,97]), max(test_matrice_sorted[,4],test_matrice_sorted[,97]), len=10000) 
p0 <- F94min(tmpx)
p1 <- F94max(tmpx)
segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "gray88", "gray88"), alpha.f=0.1))

tmpx <- seq(min(test_matrice_sorted[,14],test_matrice_sorted[,87]), max(test_matrice_sorted[,14],test_matrice_sorted[,87]), len=10000) 
p0 <- F74min(tmpx)
p1 <- F74max(tmpx)
segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "gray78", "gray78"), alpha.f=0.1))


tmpx <- seq(min(test_matrice_sorted[,26],test_matrice_sorted[,75]), max(test_matrice_sorted[,26],test_matrice_sorted[,75]), len=10000) 
p0 <- F50min(tmpx)
p1 <- F50max(tmpx)
segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "gray69", "gray68"), alpha.f=0.1))


plot(ecdf(simu_all), verticals = TRUE, do.points = FALSE, lwd=2, add=TRUE)







######



plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,1), xli=c(0,1))

plot(ecdf(as.numeric(unlist(Astayanax_CF_mutpred_pigm[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="firebrick1", lwd=3, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(Ldentata_mutpred_pigm[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="brown4", lwd=3, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(Lholguinensis_mutpred_pigm[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="orange", lwd=3, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(Astayanax_SF_mutpred_pigm[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="deeppink", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(Brotula_barbata_mutpred_pigm[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="chartreuse4", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(Lamprogrammus_mutpred_pigm[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="darkgreen", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(Carapus_acus_mutpred_pigm[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="aquamarine", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(Pygocentrus_mutpred_pigm[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="darkviolet", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(danio_rerio_mutpred_pigm[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="green", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(grahami_mutpred_pigm[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="blue4", lwd=1.4, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(anshuiensis_mutpred_pigm[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="darkred", lwd=1.4, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(rhinocerous_mutpred_pigm[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="goldenrod1", lwd=1.4, cex=0.9, pch=2)



legend("bottomright", legend=c("Simulations","Mean of 100 simulations", "Astayanax mexicanus SF (299)", "Pygocentrus nattereri (7634)","Danio rerio (9169)", "Brotula barbata (7705)", "Astyanax mexicanus CF (232)","Lucifuga gibarensis (446)", "Lucifuga dentata (834)"), col=c("gray93", "black", "deeppink", "darkviolet", "green", "chartreuse4", "firebrick1", "orange", "brown4"), lty=1, text.font=3, lwd=c(3, 3, 3, 3, 3, 3, 4, 4, 4), bty= 'n', pch=c(NA, NA, NA, NA, NA, NA, 17, 17, 17))
legend("bottomright", legend=c("Simulations","Mean of 100 simulations", "Sinocyclocheilus grahami (2639)", "Sinocyclocheilus rhinocerous (2032)","Sinocyclocheilus anshuiensis (2045)"), col=c("gray93", "black", "blue4", "goldenrod1", "darkred"), lty=1, text.font=3, lwd=c(3, 3, 3, 3, 3), bty= 'n', pch=c(NA, NA, 17, 17, 17))



dep1 <- rep("NUL", length(simu_all))

simu_all_matrix <- cbind.data.frame(dep1, dep1, simu_all)

n=0
matrix_pvalues <- matrix(nrow=8, ncol=8)
malist <- list(simu_all_matrix, Astayanax_SF_mutpred_pigm, Pygocentrus_mutpred_pigm, danio_rerio_mutpred_pigm, Brotula_barbata_mutpred_pigm, Astayanax_CF_mutpred_pigm, Lholguinensis_mutpred_pigm, Ldentata_mutpred_pigm)
tests_concat <- c()
for(i in malist){
  tests_concat <- c()
  for(j in malist){
    tests_concat <- c(tests_concat, ks.test(as.numeric(unlist(i[,3])), as.numeric(unlist(j[,3])))$p)
  }
  n=n+1
  matrix_pvalues[n,] <- tests_concat
}

cellcol[matrix_pvalues<0.05] <- "brown4"
cellcol[matrix_pvalues>0.05] <- "forestgreen"

color2D.matplot(matrix_pvalues, show.values=TRUE,cellcolors=cellcol, xlab="", ylab="",main="Kolmogorov-Smirnov Tests p-values between MutPred scores distribution", axes=FALSE, vcex=2.5, cex.main=2)



#### KERNEL DENSITIES ##

plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,5), xlim=c(0,1), cex.axis=2, cex.lab=1.5, main="Kernel densities of pigmentation genes mutpred scores", cex.main=2.3)


matrice_de_y=matrix(nrow=512)
matrice_de_x=matrix(nrow=512)
mean_distrib_simu <- c()
for(i in unique(Pigmentation_Simulations_Danio$"Simu")){
  matrice_de_x <- cbind(matrice_de_x, density(sample(as.numeric(unlist(c(Pigmentation_Simulations_Danio[Pigmentation_Simulations_Danio$Simu==i][,4]))), 1000), from=0, to=1)$x)
  matrice_de_y <- cbind(matrice_de_y, density(sample(as.numeric(unlist(c(Pigmentation_Simulations_Danio[Pigmentation_Simulations_Danio$Simu==i][,4]))), 1000), from=0, to=1)$y)
  mean_distrib_simu <- c(mean_distrib_simu, sample(as.numeric(unlist(c(Pigmentation_Simulations_Danio[Pigmentation_Simulations_Danio$Simu==i][,4]))), 1000))
}


matrice_de_x <- matrice_de_x[,-1]
matrice_de_y <- matrice_de_y[,-1]
moyenne_kernel_x <- apply(matrice_de_x, 1, mean) 
moyenne_kernel_y <- apply(matrice_de_y, 1, mean) 
max_kernel_densities_x <- apply(matrice_de_x, 1, max) 
max_kernel_density_y <- apply(matrice_de_y, 1, max)
min_kernel_density_y <- apply(matrice_de_y, 1, min)

for(i in unique(Pigmentation_Simulations_Danio$"Simu")){
  lines(density(sample(as.numeric(unlist(c(Pigmentation_Simulations_Danio[Pigmentation_Simulations_Danio$Simu==i][,4]))), 1000), from=0, to=1), col="gray88")
}


#polygon(c(max_kernel_densities_x,rev(max_kernel_densities_x)),c(max_kernel_density_y,rev(min_kernel_density_y)), col="grey88", border=NA)



sorted_coordinates <- apply(matrice_de_y, 1, sort) 
min_values <- sorted_coordinates[11,]
max_values <- sorted_coordinates[90,]
polygon(c(max_kernel_densities_x,rev(max_kernel_densities_x)),c(max_values,rev(min_values)), col="grey70", border=NA)

lines(x=moyenne_kernel_x, y=moyenne_kernel_y, col="black", lwd=2)


lines(density(as.numeric(unlist(Astayanax_CF_mutpred_pigm[,3]))), col="firebrick1", lwd=3)
lines(density(as.numeric(unlist(Ldentata_mutpred_pigm[,3]))), col="brown4", lwd=3)
lines(density(as.numeric(unlist(Lholguinensis_mutpred_pigm[,3]))), col="orange", lwd=3)
lines(density(as.numeric(unlist(Astayanax_SF_mutpred_pigm[,3]))), col="deeppink", lwd=2)
lines(density(as.numeric(unlist(Brotula_barbata_mutpred_pigm[,3]))), col="chartreuse4", lwd=2)
lines(density(as.numeric(unlist(Lamprogrammus_mutpred_pigm[,3]))), col="darkgreen", lwd=2)
lines(density(as.numeric(unlist(Carapus_acus_mutpred_pigm[,3]))),  col="aquamarine", lwd=2)
lines(density(as.numeric(unlist(Pygocentrus_mutpred_pigm[,3]))), col="darkviolet", lwd=2)
lines(density(as.numeric(unlist(danio_rerio_mutpred_pigm[,3]))), col="green", lwd=2)
lines(density(as.numeric(unlist(grahami_mutpred_pigm[,3]))), col="blue4", lwd=2)
lines(density(as.numeric(unlist(anshuiensis_mutpred_pigm[,3]))), col="cyan", lwd=2)
lines(density(as.numeric(unlist(rhinocerous_mutpred_pigm[,3]))), col="deepskyblue", lwd=2)


legend("bottomright", legend=c("Simulations","Mean of 100 simulations", "Astayanax mexicanus SF (299)", "Pygocentrus nattereri (7634)","Danio rerio (9169)", "Brotula barbata (7705)", "Astyanax mexicanus CF (232)","Lucifuga gibarensis (446)", "Lucifuga dentata (834)"), col=c("gray93", "black", "deeppink", "darkviolet", "green", "chartreuse4", "firebrick1", "orange", "brown4"), lty=1, text.font=3, lwd=c(3, 3, 3, 3, 3, 3, 4, 4, 4), bty= 'n')






###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################




library("Biostrings")
data(PAM30)


PAM5_matrix <- matrix(c(7,-12,-9,-8,-12,-9,-7,-6,-13,-10,-11,-12,-10,-14,-6,-5,-5,-24,-13,-7,-8,-8,-8,-28
                        ,-12,9,-12,-21,-13,-6,-19,-15,-6,-10,-14,-4,-9,-14,-9,-8,-12,-7,-16,-13,-14,-9,-11,-28
                        ,-9,-12,9,-2,-21,-9,-7,-8,-4,-10,-12,-5,-19,-14,-11,-4,-6,-13,-9,-14,7,-8,-8,-28
                        ,-8,-21,-2,9,-25,-8,-2,-8,-9,-13,-23,-10,-21,-25,-14,-9,-10,-26,-22,-14,7,-3,-12,-28
                        ,-12,-13,-21,-25,10,-25,-25,-15,-12,-11,-25,-25,-24,-23,-13,-7,-14,-26,-9,-11,-22,-25,-16,-28
                        ,-9,-6,-9,-8,-25,9,-3,-12,-4,-13,-10,-8,-8,-24,-8,-10,-11,-23,-22,-12,-8,7,-10,-28
                        ,-7,-19,-7,-2,-25,-3,8,-9,-11,-10,-15,-9,-12,-25,-11,-9,-12,-28,-13,-12,-3,7,-10,-28
                        ,-6,-15,-8,-8,-15,-12,-9,7,-15,-21,-16,-12,-14,-14,-12,-6,-12,-25,-24,-10,-8,-10,-10,-28
                        ,-13,-6,-4,-9,-12,-4,-11,-15,10,-16,-11,-12,-21,-11,-9,-12,-13,-12,-8,-11,-6,-6,-10,-28
                        ,-10,-10,-10,-13,-11,-13,-10,-21,-16,9,-6,-11,-5,-7,-15,-13,-7,-24,-11,-2,-11,-11,-10,-28
                        ,-11,-14,-12,-23,-25,-10,-15,-16,-11,-6,7,-14,-4,-7,-12,-14,-12,-11,-12,-7,-14,-11,-11,-28
                        ,-12,-4,-5,-10,-25,-8,-9,-12,-12,-11,-14,7,-6,-24,-12,-9,-8,-22,-14,-15,-7,-8,-10,-28
                        ,-10,-9,-19,-21,-24,-8,-12,-14,-21,-5,-4,-6,12,-9,-14,-10,-9,-23,-22,-6,-20,-10,-11,-28
                        ,-14,-14,-14,-25,-23,-24,-25,-14,-11,-7,-7,-24,-9,9,-15,-11,-14,-9,-3,-14,-16,-24,-14,-28
                        ,-6,-9,-11,-14,-13,-8,-11,-12,-9,-15,-12,-12,-14,-15,8,-6,-9,-24,-24,-11,-12,-9,-10,-28
                        ,-5,-8,-4,-9,-7,-10,-9,-6,-12,-13,-14,-9,-10,-11,-6,7,-4,-10,-12,-12,-6,-10,-8,-28
                        ,-5,-12,-6,-10,-14,-11,-12,-12,-13,-7,-12,-8,-9,-14,-9,-4,8,-23,-11,-7,-8,-11,-8,-28
                        ,-24,-7,-13,-26,-26,-23,-28,-25,-12,-24,-11,-22,-23,-9,-24,-10,-23,13,-10,-26,-15,-25,-19,-28
                        ,-13,-16,-9,-22,-9,-22,-13,-24,-8,-11,-12,-14,-22,-3,-24,-12,-11,-10,10,-12,-11,-15,-14,-28
                        ,-7,-13,-14,-14,-11,-12,-12,-10,-11,-2,-7,-15,-6,-14,-11,-12,-7,-26,-12,8,-14,-12,-9,-28
                        ,-8,-14,7,7,-22,-8,-3,-8,-6,-11,-14,-7,-20,-16,-12,-6,-8,-15,-11,-14,7,-5,-10,-28
                        ,-8,-9,-8,-3,-25,7,7,-10,-6,-11,-11,-8,-10,-24,-9,-10,-11,-25,-15,-12,-5,7,-10,-28
                        ,-8,-11,-8,-12,-16,-10,-10,-10,-10,-10,-11,-10,-11,-14,-10,-8,-8,-19,-14,-9,-10,-10,-10,-28
                        ,-28,-28,-28,-28,-28,-28,-28,-28,-28,-28,-28,-28,-28,-28,-28,-28,-28,-28,-28,-28,-28,-28,-28,1), 
                      nrow=24, byrow=TRUE)
rownames(PAM5_matrix) <- c("A" ,"R" ,"N" ,"D" ,"C" ,"Q" ,"E" ,"G" ,"H" ,"I" ,"L" ,"K" ,"M" ,"F" ,"P" ,"S" ,"T" ,"W" ,"Y" ,"V" ,"B" ,"Z" ,"X" ,"*")
colnames(PAM5_matrix) <- rownames(PAM5_matrix)




plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "- (PAM5 Score)", ylim=c(0,1), xlim=c(-13,28))

simu_all_pam5 <- c()
for(i in unique(Vision_Simulations_Danio$"Simu")){
  aa1_list <- c()
  aa2_list <- c()
  for(j in unlist(Vision_Simulations_Danio[Vision_Simulations_Danio$Simu==i][,3])){
    aa1_list <- c(aa1_list, substr(j, 1, 1))
    aa2_list <- c(aa2_list, substrRight(j, 1))
  }
  PAM5_matrix_scores <- c()
  for(x in 1:length(aa1_list)){
    PAM5_matrix_scores <- c(PAM5_matrix_scores, PAM5_matrix[aa1_list[x], aa2_list[x]])
  }
  #plot(ecdf(-PAM5_matrix_scores), verticals = TRUE, do.points = TRUE, add=TRUE, col="grey", lwd=1.4, cex=0.4, pch=19)
  #lines(density(-PAM5_matrix_scores, from=-10, to=32), col="gray", lwd=2, add=TRUE, ylim=c(0,0.20))
  simu_all_pam5 <- c(simu_all_pam5, PAM5_matrix_scores)
}




plot(density(-PAM5_matrix_scores, from=-10, to=32), col="gray", lwd=2, add=TRUE, ylim=c(0,0.20))


malist <- list(-simu_all_pam5, -PAM5_matrix_scores)

multhist(malist, freq=FALSE, col=c("gray", "forestgreen"), ylim=c(0, 0.2), xlab="-PAM5 scores", ylab="Frequency")


###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(Astayanax_CF_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores <- c(PAM5_matrix_scores, PAM5_matrix[aa1_list[i], aa2_list[i]])
}

plot(ecdf(PAM5_matrix_scores), verticals = TRUE, do.points = TRUE, add=TRUE, col="firebrick1", lwd=1.4, cex=0.4, pch=19)


###################################################################################

aa1_list <- c()
aa2_list <- c()

for(i in unlist(Astayanax_SF_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores <- c(PAM5_matrix_scores, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(PAM5_matrix_scores), verticals = TRUE, do.points = TRUE, add=TRUE, col="deeppink", lwd=1.4, cex=0.4, pch=19)


###################################################################################

aa1_list <- c()
aa2_list <- c()

for(i in unlist(Ldentata_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores <- c(PAM5_matrix_scores, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(PAM5_matrix_scores), verticals = TRUE, do.points = TRUE, add=TRUE, col="brown4", lwd=1.4, cex=0.4, pch=19)


###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(Lholguinensis_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores <- c(PAM5_matrix_scores, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(PAM5_matrix_scores), verticals = TRUE, do.points = TRUE, add=TRUE, col="orange", lwd=1.4, cex=0.4, pch=19)


###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(Brotula_barbata_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores <- c(PAM5_matrix_scores, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(PAM5_matrix_scores), verticals = TRUE, do.points = TRUE, add=TRUE, col="green", lwd=1.4, cex=0.4, pch=19)

###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(Lamprogrammus_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores <- c(PAM5_matrix_scores, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(PAM5_matrix_scores), verticals = TRUE, do.points = TRUE, add=TRUE, col="chartreuse4", lwd=1.4, cex=0.4, pch=19)


###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(Carapus_acus_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores <- c(PAM5_matrix_scores, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(PAM5_matrix_scores), verticals = TRUE, do.points = TRUE, add=TRUE, col="darkgreen", lwd=1.4, cex=0.4, pch=19)

###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(Pygocentrus_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores <- c(PAM5_matrix_scores, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(PAM5_matrix_scores), verticals = TRUE, do.points = TRUE, add=TRUE, col="aquamarine", lwd=1.4, cex=0.4, pch=19)

###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(danio_rerio_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores <- c(PAM5_matrix_scores, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(-PAM5_matrix_scores), verticals = TRUE, do.points = TRUE, add=TRUE, col="forestgreen", lwd=1.4, cex=0.4, pch=19)


lines(density(-PAM5_matrix_scores, from=-10, to=32), col="forestgreen", lwd=2, add=TRUE)


###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(grahami_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores <- c(PAM5_matrix_scores, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(PAM5_matrix_scores), verticals = TRUE, do.points = TRUE, add=TRUE, col="forestgreen", lwd=1.4, cex=0.4, pch=19)


###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(anshuiensis_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores <- c(PAM5_matrix_scores, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(PAM5_matrix_scores), verticals = TRUE, do.points = TRUE, add=TRUE, col="blue4", lwd=1.4, cex=0.4, pch=19)


###################################################################################


aa1_list <- c()
aa2_list <- c()

for(i in unlist(rhinocerous_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores <- c(PAM5_matrix_scores, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


plot(ecdf(PAM5_matrix_scores), verticals = TRUE, do.points = TRUE, add=TRUE, col="deepskyblue", lwd=1.4, cex=0.4, pch=19)





#### Correlation scores grantham, pam5 et mutpred ####
### Les deux Astyanax ###



library(corrplot)
library(ggpubr)



aa1_list <- c()
aa2_list <- c()

for(i in unlist(Astayanax_CF_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores_CF <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores_CF <- c(PAM5_matrix_scores_CF, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


grantham_list_CF <- c()
for(i in 1:length(aa1_list)){
  grantham_list_CF <- c(grantham_list_CF, Grantham_matrix[aa1_list[i], aa2_list[i]])
}


aa1_list <- c()
aa2_list <- c()

for(i in unlist(Astayanax_SF_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores_SF <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores_SF <- c(PAM5_matrix_scores_SF, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


grantham_list_SF <- c()
for(i in 1:length(aa1_list)){
  grantham_list_SF <- c(grantham_list_SF, Grantham_matrix[aa1_list[i], aa2_list[i]])
}



aa1_list <- c()
aa2_list <- c()

for(i in unlist(Pygocentrus_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores_Pygo <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores_Pygo <- c(PAM5_matrix_scores_Pygo, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


grantham_list_Pygo <- c()
for(i in 1:length(aa1_list)){
  grantham_list_Pygo <- c(grantham_list_Pygo, Grantham_matrix[aa1_list[i], aa2_list[i]])
}








Cave_rep <- rep("Astyanax mexicanus Cave", 54)
Surface_rep <- rep("Astyanax mexicanus Surface", 45)
Pygo_rep <- rep("Pygocentrus nattereri", length(as.numeric(unlist(Pygocentrus_mutpred[,3]))))

Mutpred_list_CF <- as.numeric(unlist(Astayanax_CF_mutpred[,3]))
Mutpred_list_SF <- as.numeric(unlist(Astayanax_SF_mutpred[,3]))
MutPred_list_Pygo <- as.numeric(unlist(Pygocentrus_mutpred[,3]))


matable_cor_CF <- cbind.data.frame(Mutpred_list_CF, PAM5_matrix_scores_CF, grantham_list_CF, Cave_rep)
colnames(matable_cor_CF) <- c("MutPred", "PAM5", "Grantham", "Species")
matable_cor_SF <- cbind.data.frame(Mutpred_list_SF, PAM5_matrix_scores_SF, grantham_list_SF, Surface_rep)
colnames(matable_cor_SF) <- c("MutPred", "PAM5", "Grantham", "Species")
matable_cor_Pygo <- cbind.data.frame(MutPred_list_Pygo, PAM5_matrix_scores_Pygo, grantham_list_Pygo, Pygo_rep)
colnames(matable_cor_Pygo) <- c("MutPred", "PAM5", "Grantham", "Species")

#matable_cor <- rbind(matable_cor_CF, matable_cor_SF, matable_cor_Pygo)
matable_cor <- rbind(matable_cor_CF, matable_cor_SF)


cor.test(Mutpred_list_CF, grantham_list_CF, method="kendall")


ggscatter(matable_cor, x="Grantham", y="PAM5", color="Species", add="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="pearson", xlab="Grantham Score", ylab="PAM5 Score")
ggscatter(matable_cor, x="MutPred", y="PAM5", add="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="pearson", xlab="MutPred Score", ylab="PAM5 Score")
ggscatter(matable_cor, x="Grantham", y="MutPred", color="Species", add="reg.line", cor.method="kendall",  cor.coef=TRUE, xlab="Grantham Score", ylab="MutPred Score")

sp <- ggscatter(matable_cor, x = "Grantham", y = "MutPred",
                cor.method="kendall",
                add = "reg.line",
                cor.coeff.args = list(method = "kendall", label.x = 3, label.sep = "\n"), 
                conf.int = TRUE,                
                color = "Species", palette = "jco", # Color by groups "cyl"
                shape = "Species"                   # Change point shape by groups "cyl"
)+
  stat_cor(aes(color = Species),method="kendall", label.x = 3)       # Add correlation coefficient
sp 

xplot <- ggdensity(matable_cor, "Grantham", fill = "Species",
                   palette = "jco")
yplot <- ggdensity(matable_cor, "MutPred", fill = "Species", 
                   palette = "jco")+
  rotate()

yplot <- yplot + clean_theme() 
xplot <- xplot + clean_theme()

ggarrange(xplot, NULL, sp, yplot, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)




wilcox.test(Mutpred_list_CF, MutPred_list_Pygo)
wilcox.test(Mutpred_list_SF, MutPred_list_Pygo)
wilcox.test(Mutpred_list_CF, Mutpred_list_SF)



## Les deux Lucifuga ##



aa1_list <- c()
aa2_list <- c()

for(i in unlist(Ldentata_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores_Ldent <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores_Ldent <- c(PAM5_matrix_scores_Ldent, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


grantham_list_Ldent <- c()
for(i in 1:length(aa1_list)){
  grantham_list_Ldent <- c(grantham_list_Ldent, Grantham_matrix[aa1_list[i], aa2_list[i]])
}


aa1_list <- c()
aa2_list <- c()

for(i in unlist(Lholguinensis_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores_Lholg <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores_Lholg <- c(PAM5_matrix_scores_Lholg, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


grantham_list_Lholg <- c()
for(i in 1:length(aa1_list)){
  grantham_list_Lholg <- c(grantham_list_Lholg, Grantham_matrix[aa1_list[i], aa2_list[i]])
}


##
aa1_list <- c()
aa2_list <- c()

for(i in unlist(Brotula_barbata_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores_Bbarbata <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores_Bbarbata <- c(PAM5_matrix_scores_Bbarbata, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


grantham_list_Bbarbata <- c()
for(i in 1:length(aa1_list)){
  grantham_list_Bbarbata <- c(grantham_list_Bbarbata, Grantham_matrix[aa1_list[i], aa2_list[i]])
}
##

aa1_list <- c()
aa2_list <- c()

for(i in unlist(Carapus_acus_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores_Cacus <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores_Cacus <- c(PAM5_matrix_scores_Cacus, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


grantham_list_Cacus <- c()
for(i in 1:length(aa1_list)){
  grantham_list_Cacus <- c(grantham_list_Cacus, Grantham_matrix[aa1_list[i], aa2_list[i]])
}

##

aa1_list <- c()
aa2_list <- c()

for(i in unlist(Lamprogrammus_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores_Lexutus <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores_Lexutus <- c(PAM5_matrix_scores_Lexutus, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


grantham_list_Lexutus <- c()
for(i in 1:length(aa1_list)){
  grantham_list_Lexutus <- c(grantham_list_Lexutus, Grantham_matrix[aa1_list[i], aa2_list[i]])
}

##


Lexutus_rep <- rep("Lamprogrammus exutus", length(unlist(Lamprogrammus_mutpred[,2])))
Cacus_rep <- rep("Carapus acus", length(unlist(Carapus_acus_mutpred[,2])))
Brotula_rep <- rep("Brotula barbata", length(unlist(Brotula_barbata_mutpred[,2])))
Cave_rep <- rep("Lucifuga dentata", length(unlist(Ldentata_mutpred[,2])))
Surface_rep <- rep("Lucifuga holguinensis", length(unlist(Lholguinensis_mutpred[,2])))




Mutpred_list_Ldent <- as.numeric(unlist(Ldentata_mutpred[,3]))
Mutpred_list_Lholg <- as.numeric(unlist(Lholguinensis_mutpred[,3]))
Mutpred_list_Bbarbata <- as.numeric(unlist(Brotula_barbata_mutpred[,3]))
Mutpred_list_Cacus <- as.numeric(unlist(Carapus_acus_mutpred[,3]))
Mutpred_list_Lexutus <- as.numeric(unlist(Lamprogrammus_mutpred[,3]))


matable_cor_Ldent <- cbind.data.frame(Mutpred_list_Ldent, PAM5_matrix_scores_Ldent, grantham_list_Ldent, Cave_rep)
colnames(matable_cor_Ldent) <- c("MutPred", "PAM5", "Grantham", "Species")
matable_cor_Lholg <- cbind.data.frame(Mutpred_list_Lholg, PAM5_matrix_scores_Lholg, grantham_list_Lholg, Surface_rep)
colnames(matable_cor_Lholg) <- c("MutPred", "PAM5", "Grantham", "Species")
matable_cor_Bbarbata <- cbind.data.frame(Mutpred_list_Bbarbata, PAM5_matrix_scores_Bbarbata, grantham_list_Bbarbata, Brotula_rep)
colnames(matable_cor_Bbarbata) <- c("MutPred", "PAM5", "Grantham", "Species")
matable_cor_Cacus <- cbind.data.frame(Mutpred_list_Cacus, PAM5_matrix_scores_Cacus, grantham_list_Cacus, Cacus_rep)
colnames(matable_cor_Cacus) <- c("MutPred", "PAM5", "Grantham", "Species")
matable_cor_Lexutus <- cbind.data.frame(Mutpred_list_Lexutus, PAM5_matrix_scores_Lexutus, grantham_list_Lexutus, Lexutus_rep)
colnames(matable_cor_Lexutus) <- c("MutPred", "PAM5", "Grantham", "Species")
matable_cor <- rbind(matable_cor_Ldent, matable_cor_Lholg, matable_cor_Bbarbata, matable_cor_Cacus, matable_cor_Lexutus)
matable_cor <- rbind(matable_cor_Ldent, matable_cor_Lholg)


cor.test(Mutpred_list_Ldent, grantham_list_Ldent, method="kendall")


ggscatter(matable_cor, x="Grantham", y="PAM5", color="Species", add="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="pearson", xlab="Grantham Score", ylab="PAM5 Score")
ggscatter(matable_cor, x="MutPred", y="PAM5", add="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="pearson", xlab="MutPred Score", ylab="PAM5 Score")
ggscatter(matable_cor, x="Grantham", y="MutPred", color="Species", add="reg.line", cor.method="kendall",  cor.coef=TRUE, xlab="Grantham Score", ylab="MutPred Score")

sp <- ggscatter(matable_cor, x = "Grantham", y = "MutPred",
                cor.method="kendall",
                #size=0.6,
                add = "reg.line",
                cor.coeff.args = list(method = "kendall", label.x = 3, label.sep = "\n"), 
                conf.int = TRUE,                
                color = "Species", palette = "jco", # Color by groups "cyl"
                shape = "Species"                   # Change point shape by groups "cyl"
)+
  stat_cor(aes(color = Species),method="kendall", label.x = 3)       # Add correlation coefficient
sp 

xplot <- ggdensity(matable_cor, "Grantham", fill = "Species",
                   palette = "jco")
yplot <- ggdensity(matable_cor, "MutPred", fill = "Species", 
                   palette = "jco")+
  rotate()

yplot <- yplot + clean_theme() 
xplot <- xplot + clean_theme()

ggarrange(xplot, NULL, sp, yplot, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)


wilcox.test(Mutpred_list_Ldent, Mutpred_list_Lholg)
wilcox.test(Mutpred_list_Ldent, Mutpred_list_Bbarbata)
wilcox.test(Mutpred_list_Lholg, Mutpred_list_Bbarbata)
wilcox.test(Mutpred_list_Cacus, Mutpred_list_Bbarbata)
wilcox.test(Mutpred_list_Lexutus, Mutpred_list_Bbarbata)
wilcox.test(Mutpred_list_Lexutus, Mutpred_list_Cacus)




##################################################################
##################################################################
##################################################################
##################################################################
##################################################################



aa1_list <- c()
aa2_list <- c()

for(i in unlist(anshuiensis_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores_ansh <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores_ansh <- c(PAM5_matrix_scores_ansh, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


grantham_list_ansh <- c()
for(i in 1:length(aa1_list)){
  grantham_list_ansh <- c(grantham_list_ansh, Grantham_matrix[aa1_list[i], aa2_list[i]])
}


ansh_rep <- rep("Sinocyclocheilus anshuiensis", length(unlist(anshuiensis_mutpred[,2])))
Mutpred_list_ansh <- as.numeric(unlist(anshuiensis_mutpred[,3]))


##

aa1_list <- c()
aa2_list <- c()

for(i in unlist(grahami_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores_ansh <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores_ansh <- c(PAM5_matrix_scores_ansh, PAM5_matrix[aa1_list[i], aa2_list[i]])
}


grantham_list_ansh <- c()
for(i in 1:length(aa1_list)){
  grantham_list_ansh <- c(grantham_list_ansh, Grantham_matrix[aa1_list[i], aa2_list[i]])
}


ansh_rep <- rep("Sinocyclocheilus grahami", length(unlist(grahami_mutpred[,2])))
Mutpred_list_ansh <- as.numeric(unlist(grahami_mutpred[,3]))






grahami_mutpred <- list.files(path = ".", recursive = T, pattern = "grahami_vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

rhinocerous_mutpred <- list.files(path = ".", recursive = T, pattern = "rhinocerous_vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

danio_rerio_mutpred <- list.files(path = ".", recursive = T, pattern = "Danio_rerio_vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))






################################################################################


sp <- ggscatter(matable_cor, x = "Grantham", y = "PAM5",
                cor.method="kendall",
                #size=0.6,
                add = "reg.line",
                cor.coeff.args = list(method = "kendall", label.x = 3, label.sep = "\n"), 
                conf.int = TRUE,                
                color = "Species", palette = "jco", # Color by groups "cyl"
                shape = "Species"                   # Change point shape by groups "cyl"
)+
  stat_cor(aes(color = Species),method="kendall", label.x = 3)       # Add correlation coefficient
sp 

xplot <- ggdensity(matable_cor, "Grantham", fill = "Species",
                   palette = "jco")
yplot <- ggdensity(matable_cor, "PAM5", fill = "Species", 
                   palette = "jco")+
  rotate()

yplot <- yplot + clean_theme() 
xplot <- xplot + clean_theme()

ggarrange(xplot, NULL, sp, yplot, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)



################################################################################
#################### Histograms comparisons ####################################
################################################################################

simu_all_pam5 <- c()
for(i in unique(Vision_Simulations_Danio$"Simu")){
  aa1_list <- c()
  aa2_list <- c()
  for(j in unlist(Vision_Simulations_Danio[Vision_Simulations_Danio$Simu==i][,3])){
    aa1_list <- c(aa1_list, substr(j, 1, 1))
    aa2_list <- c(aa2_list, substrRight(j, 1))
  }
  PAM5_matrix_scores <- c()
  for(x in 1:length(aa1_list)){
    PAM5_matrix_scores <- c(PAM5_matrix_scores, PAM5_matrix[aa1_list[x], aa2_list[x]])
  }
  #plot(ecdf(-PAM5_matrix_scores), verticals = TRUE, do.points = TRUE, add=TRUE, col="grey", lwd=1.4, cex=0.4, pch=19)
  #lines(density(-PAM5_matrix_scores, from=-10, to=32), col="gray", lwd=2, add=TRUE, ylim=c(0,0.20))
  simu_all_pam5 <- c(simu_all_pam5, PAM5_matrix_scores)
}


simu_all_grantham <- c()
for(i in unique(Vision_Simulations_Danio$"Simu")){
  aa1_list <- c()
  aa2_list <- c()
  for(j in unlist(Vision_Simulations_Danio[Vision_Simulations_Danio$Simu==i][,3])){
    aa1_list <- c(aa1_list, substr(j, 1, 1))
    aa2_list <- c(aa2_list, substrRight(j, 1))
  }
  grantham_list <- c()
  for(x in 1:length(aa1_list)){
    grantham_list <- c(grantham_list, Grantham_matrix[aa1_list[x], aa2_list[x]])
  }
  simu_all_grantham <- c(simu_all_grantham, grantham_list)
}



aa1_list <- c()
aa2_list <- c()

for(i in unlist(danio_rerio_mutpred[,2])){
  aa1_list <- c(aa1_list, substr(i, 1, 1))
  aa2_list <- c(aa2_list, substrRight(i, 1))
  
}

PAM5_matrix_scores <- c()
grantham_list <- c()
for(i in 1:length(aa1_list)){
  PAM5_matrix_scores <- c(PAM5_matrix_scores, PAM5_matrix[aa1_list[i], aa2_list[i]])
  grantham_list <- c(grantham_list, Grantham_matrix[aa1_list[i], aa2_list[i]])
  
}


Danio_rerio_mutpred_scores <- as.numeric(unlist(danio_rerio_mutpred[,3]))

cat1_danio <- sum(Danio_rerio_mutpred_scores > 0 & Danio_rerio_mutpred_scores <= 0.1)
cat2_danio <- sum(Danio_rerio_mutpred_scores > 0.1 & Danio_rerio_mutpred_scores <= 0.2)
cat3_danio <- sum(Danio_rerio_mutpred_scores > 0.2 & Danio_rerio_mutpred_scores <= 0.3)
cat4_danio <- sum(Danio_rerio_mutpred_scores > 0.3 & Danio_rerio_mutpred_scores <= 0.4)
cat5_danio <- sum(Danio_rerio_mutpred_scores > 0.4 & Danio_rerio_mutpred_scores <= 0.5)
cat6_danio <- sum(Danio_rerio_mutpred_scores > 0.5 & Danio_rerio_mutpred_scores <= 0.6)
cat7_danio <- sum(Danio_rerio_mutpred_scores > 0.6 & Danio_rerio_mutpred_scores <= 0.7)
cat8_danio <- sum(Danio_rerio_mutpred_scores > 0.7 & Danio_rerio_mutpred_scores <= 0.8)
cat9_danio <- sum(Danio_rerio_mutpred_scores > 0.8 & Danio_rerio_mutpred_scores <= 0.9)
cat10_danio <- sum(Danio_rerio_mutpred_scores > 0.9 & Danio_rerio_mutpred_scores <= 1)

cat1_neutral <- sum(simu_all > 0 & simu_all <= 0.1)
cat2_neutral <- sum(simu_all > 0.1 & simu_all <= 0.2)
cat3_neutral <- sum(simu_all > 0.2 & simu_all <= 0.3)
cat4_neutral <- sum(simu_all > 0.3 & simu_all <= 0.4)
cat5_neutral <- sum(simu_all > 0.4 & simu_all <= 0.5)
cat6_neutral <- sum(simu_all > 0.5 & simu_all <= 0.6)
cat7_neutral <- sum(simu_all > 0.6 & simu_all <= 0.7)
cat8_neutral <- sum(simu_all > 0.7 & simu_all <= 0.8)
cat9_neutral <- sum(simu_all > 0.8 & simu_all <= 0.9)
cat10_neutral <- sum(simu_all > 0.9 & simu_all <= 1)


neutral_mutpred_hist <- c(cat1_neutral,cat2_neutral,cat3_neutral,cat4_neutral,cat5_neutral,cat6_neutral,cat7_neutral,cat8_neutral,cat9_neutral,cat10_neutral)
danio_mutpred_hist <- c(cat1_danio,cat2_danio,cat3_danio,cat4_danio,cat5_danio,cat6_danio,cat7_danio,cat8_danio,cat9_danio,cat10_danio)

neutral_mutpred_hist <- neutral_mutpred_hist/sum(neutral_mutpred_hist)
danio_mutpred_hist <- danio_mutpred_hist/sum(danio_mutpred_hist)
names_cat <- c("0.1", "0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1")

par(mfrow=c(1,3), xpd=TRUE)

malist_PAM5 <- list(-simu_all_pam5, -PAM5_matrix_scores)
malist_grantham <- list(simu_all_grantham, grantham_list)
malist_mutpred <- rbind(neutral_mutpred_hist, danio_mutpred_hist)

multhist(malist_grantham, freq=FALSE, col=c("gray", "forestgreen"), ylim=c(0, 0.03), xlab="Grantham scores", ylab="Frequency")
multhist(malist_PAM5, freq=FALSE, col=c("gray", "forestgreen"), ylim=c(0, 0.2), xlab="-PAM5 scores", ylab="Frequency")
barplot(malist_mutpred, names.arg=names_cat, xlab="MutPred scores", ylab="Frequency", beside=T, col=c("gray", "forestgreen"), ylim=c(0, 0.3))

legend(9.6,0.3,c("Neutral mutations", "Danio rerio mutations"), pch = 15, bty="n", col=c("gray", "forestgreen"), cex=1.5)



#### Test bootstrap sur Danio vision ####


test_matrice_bootstrap <- matrix(nrow=54)
for(i in 1:10000){
  test_matrice_bootstrap <- cbind(test_matrice_bootstrap, sample(Danio_rerio_mutpred_scores, 54))
}
test_matrice_bootstrap <- test_matrice_bootstrap[,-1]


plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,1), xli=c(0,1), main="Empirical distribution of MutPred scores for vision genes mutations")


for(i in 1:10000){
  plot(ecdf(test_matrice_bootstrap[,i]), verticals = TRUE, do.points = FALSE, lwd=1, add=TRUE, col="darkseagreen1")
}

plot(ecdf(Danio_rerio_mutpred_scores), verticals = TRUE, do.points = FALSE, lwd=2, add=TRUE, col="green")




## COMPARAISON 3 DATASETS SIMULATIONS ##

plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,1), xli=c(0,1), main="Empirical distribution of MutPred scores for simulated mutations")

test_matrice <- matrix(nrow=300)
for(i in unique(Vision_Simulations_Danio$"Simu")){
  test_matrice <- cbind(test_matrice, sample(as.numeric(unlist(c(Vision_Simulations_Danio[Vision_Simulations_Danio$Simu==i][,4]))), 300))
}
test_matrice <- test_matrice[,-1]

simu_all <- c()
for(i in 1:100){
  simu_all <- c(simu_all, as.numeric(test_matrice[,i]))
}

for(i in 1:100){
  plot(ecdf(test_matrice[,i]), verticals = TRUE, do.points = FALSE, lwd=1, add=TRUE, col="lightcoral")
}

plot(ecdf(simu_all), verticals = TRUE, do.points = FALSE, lwd=2, add=TRUE, col="red3")




test_matrice <- matrix(nrow=200)
for(i in unique(Circadian_Simulations_Danio$"Simu")){
  test_matrice <- cbind(test_matrice, sample(as.numeric(unlist(c(Circadian_Simulations_Danio[Circadian_Simulations_Danio$Simu==i][,4]))), 200))
}
test_matrice <- test_matrice[,-1]

simu_all <- c()
for(i in 1:100){
  simu_all <- c(simu_all, as.numeric(test_matrice[,i]))
}

for(i in 1:100){
  plot(ecdf(test_matrice[,i]), verticals = TRUE, do.points = FALSE, lwd=1, add=TRUE, col="darkseagreen1")
}

plot(ecdf(simu_all), verticals = TRUE, do.points = FALSE, lwd=2, add=TRUE, col="forestgreen")

test_matrice <- matrix(nrow=1000)
for(i in unique(Pigmentation_Simulations_Danio$"Simu")){
  test_matrice <- cbind(test_matrice, sample(as.numeric(unlist(c(Pigmentation_Simulations_Danio[Pigmentation_Simulations_Danio$Simu==i][,4]))), 1000))
}
test_matrice <- test_matrice[,-1]

simu_all <- c()
for(i in 1:100){
  simu_all <- c(simu_all, as.numeric(test_matrice[,i]))
}

for(i in 1:100){
  plot(ecdf(test_matrice[,i]), verticals = TRUE, do.points = FALSE, lwd=1, add=TRUE, col="skyblue1")
}

plot(ecdf(simu_all), verticals = TRUE, do.points = FALSE, lwd=2, add=TRUE, col="royalblue4")

legend("bottomright", legend=c("Vision", "Pigmentation", "Circadian clock"), title="100 Simulations", lty=1, col=c("lightcoral","skyblue1" ,"darkseagreen1"), lwd=2, bty= 'n')



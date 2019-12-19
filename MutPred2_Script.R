library("boot")
library(plotrix)
library(data.table)
library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(data.table)

read_plus <- function(flnm) {
  fread(file = flnm, fill = T, sep="\t") 
}


# Import Mutpred2 results for vision, circadian and pigmentation mutations #

setwd("~/Desktop/Maximum_Likelihood_Mutpred_Pipeline/Script_Mutpred2/")


Astayanax_CF_mutpred <- list.files(path = ".", recursive = T, pattern = "Astyanax_CF_Vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))


Astayanax_SF_mutpred <- list.files(path = ".", recursive = T, pattern = "Astyanax_Surface_vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))

Brotula_barbata_mutpred <- list.files(path = ".", recursive = T, pattern = "Brotula_barbata_vision.results.parsed_scores", full.names = F) %>% # list files
  map_df(~read_plus(.))


Pygocentrus_mutpred <- list.files(path = ".", recursive = T, pattern = "Pygocentrus_nattereri_vision.results.parsed_scores", full.names = F) %>% # list files
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


Vision_Simulations_Danio <- list.files(path = ".", recursive = T, pattern = "Danio_rerio_Simulations_Vision.tsv", full.names = F) %>% # list files
  map_df(~read_plus(.))

Vision_Simulations_Pygocentrus <- list.files(path = ".", recursive = T, pattern = "Pygocentrus_Simulations_Vision.tsv", full.names = F) %>% # list files
  map_df(~read_plus(.))

Vision_Simulations_Brotula <- list.files(path = ".", recursive = T, pattern = "Brotula_barbata_Simulations_Vision.tsv", full.names = F) %>% # list files
  map_df(~read_plus(.))


colnames(Vision_Simulations_Danio) <- c("Simu", "gene", "Substitution", "MutPred Score")
colnames(Vision_Simulations_Pygocentrus) <- c("Simu", "gene", "Substitution", "MutPred Score")
colnames(Vision_Simulations_Brotula) <- c("Simu", "gene", "Substitution", "MutPred Score")



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

Circadian_Simulations_Danio <- list.files(path = ".", recursive = T, pattern = "Danio_rerio_Simulations_Circadian.tsv", full.names = F) %>% # list files
  map_df(~read_plus(.))

colnames(Circadian_Simulations_Danio) <- c("Simu", "gene", "Substitution", "MutPred Score")


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


Pigmentation_Simulations_Danio <- list.files(path = ".", recursive = T, pattern = "Danio_rerio_Simulations_Pigmentation.tsv", full.names = F) %>% # list files
  map_df(~read_plus(.))

colnames(Pigmentation_Simulations_Danio) <- c("Simu", "gene", "Substitution", "MutPred Score")


# Example of how to plot an empirical distributions (ECDF) #

test_matrice <- matrix(nrow=300)
for(i in unique(Vision_Simulations_Danio$"Simu")){
  test_matrice <- cbind(test_matrice, sample(as.numeric(unlist(c(Vision_Simulations_Danio[Vision_Simulations_Danio$Simu==i][,4]))), 300))
}
test_matrice <- test_matrice[,-1]

simu_all <- c()
for(i in 1:100){
  simu_all <- c(simu_all, as.numeric(test_matrice[,i]))
}

plot(0, bty = 'n', pch = '', ylab = "Cumulative frequency", xlab = "MutPred Score", ylim=c(0,1), xlim=c(0,1), main="Eye genes", cex.lab=1.5, cex.axis=1.7, cex.main=2.3)

for(i in 1:100){
  plot(ecdf(test_matrice[,i]), verticals = TRUE, do.points = FALSE, lwd=1, add=TRUE, col="gray93")
}

plot(ecdf(simu_all), verticals = TRUE, do.points = FALSE, lwd=2, add=TRUE)

plot(ecdf(as.numeric(unlist(Astayanax_CF_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="firebrick1", lwd=3, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(Ldentata_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="brown4", lwd=3, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(Lholguinensis_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="orange", lwd=3, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(Astayanax_SF_mutpred[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="deeppink", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(Brotula_barbata_mutpred[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="chartreuse4", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(Pygocentrus_mutpred[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="darkviolet", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(danio_rerio_mutpred[,3]))), verticals = TRUE, do.points = FALSE, add=TRUE, col="green", lwd=1.4, cex=0.4, pch=19)
plot(ecdf(as.numeric(unlist(grahami_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="blue4", lwd=0.8, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(anshuiensis_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="darkred", lwd=0.8, cex=0.9, pch=2)
plot(ecdf(as.numeric(unlist(rhinocerous_mutpred[,3]))), verticals = TRUE, do.points = TRUE, add=TRUE, col="goldenrod1", lwd=0.8, cex=0.9, pch=2)


# Example of how to plot a kernel density distribution #


plot(0, bty = 'n', pch = '', ylab = "Density", xlab = "MutPred Score", ylim=c(0,3.2), xlim=c(0,1), cex.axis=2, cex.lab=1.5, main="Kernel densities of vision genes mutpred scores", cex.main=2.3)


matrice_de_y=matrix(nrow=512)
matrice_de_x=matrix(nrow=512)
mean_distrib_simu <- c()
for(i in unique(Vision_Simulations_Danio$"Simu")){
  matrice_de_x <- cbind(matrice_de_x, density(sample(as.numeric(unlist(c(Vision_Simulations_Danio[Vision_Simulations_Danio$Simu==i][,4]))), 300), from=0, to=1)$x)
  matrice_de_y <- cbind(matrice_de_y, density(sample(as.numeric(unlist(c(Vision_Simulations_Danio[Vision_Simulations_Danio$Simu==i][,4]))), 300), from=0, to=1)$y)
  mean_distrib_simu <- c(mean_distrib_simu, sample(as.numeric(unlist(c(Vision_Simulations_Danio[Vision_Simulations_Danio$Simu==i][,4]))), 300))
}



for(i in unique(Vision_Simulations_Danio$"Simu")){
  lines(density(sample(as.numeric(unlist(c(Vision_Simulations_Danio[Vision_Simulations_Danio$Simu==i][,4]))), 300), from=0, to=1), col="gray88")
}

matrice_de_x <- matrice_de_x[,-1]
matrice_de_y <- matrice_de_y[,-1]
moyenne_kernel_x <- apply(matrice_de_x, 1, mean) 
moyenne_kernel_y <- apply(matrice_de_y, 1, mean) 
lines(x=moyenne_kernel_x, y=moyenne_kernel_y, col="black", lwd=2)


lines(density(as.numeric(unlist(Astayanax_CF_mutpred[,3])), from=0, to=1), col="firebrick1", lwd=3)
lines(density(as.numeric(unlist(Ldentata_mutpred[,3])), from=0, to=1), col="brown4", lwd=3)
lines(density(as.numeric(unlist(Lholguinensis_mutpred[,3])), from=0, to=1), col="orange", lwd=3)
lines(density(as.numeric(unlist(Astayanax_SF_mutpred[,3])), from=0, to=1), col="deeppink", lwd=2)
lines(density(as.numeric(unlist(Brotula_barbata_mutpred[,3])), from=0, to=1), col="chartreuse4", lwd=2)
lines(density(as.numeric(unlist(Pygocentrus_mutpred[,3])), from=0, to=1), col="darkviolet", lwd=2)
lines(density(as.numeric(unlist(danio_rerio_mutpred[,3])), from=0, to=1), col="green", lwd=2)
lines(density(as.numeric(unlist(grahami_mutpred[,3])), from=0, to=1), col="blue4", lwd=2)
lines(density(as.numeric(unlist(anshuiensis_mutpred[,3])), from=0, to=1), col="darkred", lwd=2)
lines(density(as.numeric(unlist(rhinocerous_mutpred[,3])), from=0, to=1), col="goldenrod1", lwd=2)


# Example of how to compute k.s tests between distributions #


library(RColorBrewer)

dep1 <- rep("NUL", length(simu_all))

simu_all_matrix <- cbind.data.frame(dep1, dep1, simu_all)

n=0
matrix_pvalues <- matrix(nrow=5, ncol=5)
malist <- list(simu_all_matrix, danio_rerio_mutpred, grahami_mutpred, rhinocerous_mutpred, anshuiensis_mutpred)
tests_concat <- c()
for(i in malist){
  tests_concat <- c()
  for(j in malist){
    tests_concat <- c(tests_concat, ks.test(as.numeric(unlist(i[,3])), as.numeric(unlist(j[,3])))$p)
  }
  n=n+1
  matrix_pvalues[n,] <- tests_concat
}


cellcol<-matrix(rep("#000000",25),nrow=5)
cellcol[matrix_pvalues<0.05] <- "brown4"
cellcol[matrix_pvalues>0.05] <- "forestgreen"

color2D.matplot(matrix_pvalues, show.values=TRUE,cellcolors=cellcol, xlab="", ylab="",main="Kolmogorov-Smirnov Tests p-values between MutPred scores distribution", axes=FALSE, vcex=2.5, cex.main=2)


# Example of how to compute grantham and pam scores for each mutations #


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

#Compute grantham scores for simulated mutations

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
  plot(ecdf(grantham_list), verticals = TRUE, do.points = TRUE, add=TRUE, col="grey", lwd=1.4, cex=0.4, pch=19)
  simu_all_grantham <- c(simu_all_grantham, grantham_list)
}


#Compute grantham scores for observed mutations (Danio rerio) and plot


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


#Same thing with PAM5

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

#Compute PAM5 scores for simulated mutations


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
  plot(ecdf(-PAM5_matrix_scores), verticals = TRUE, do.points = TRUE, add=TRUE, col="grey", lwd=1.4, cex=0.4, pch=19)
  simu_all_pam5 <- c(simu_all_pam5, PAM5_matrix_scores)
  }


#Compute PAM5 scores for observed mutations

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




### Example of an admixture simulation ###

Neutral_mutations <- simu_all #Simultated mutations represent mutations that arise under neutral evolution
#Mutations observed in surface fishes represent mutations arising under selective constraint
Positiv_mutations <- as.numeric(unlist(danio_rerio_mutpred[,3]))

Ldentata_distrib <- as.numeric(unlist(Ldentata_mutpred[,3]))
Amex_CF_distrib <- as.numeric(unlist(Astayanax_CF_mutpred[,3]))
Lholguin_distrib <- as.numeric(unlist(Lholguinensis_mutpred[,3]))


#We make 1000 simulations on which we mix those two type of mutations (for a total of 1000 mutations)

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







# Corresponding authors :Didier Casane (didier.casane@egce.cnrs-gif.fr) and Maxime Policarpo (maxime.policarpo@egce.cnrs-gif.fr)
# FIRST PART : SUPPLEMENTAL 1
# SECOND PART : SUPPLEMENTAL 2

############################################################## Stops ################################################################################ 


#Stop positions relative to the associated gene length (full length = 100%)
stop_position=c(45.4,48,95,16.57,39.19,51.74,80.2,79.1,28.64,10,36.72,2.08,18.4,79.26,2.97,27.71,43.6,30.9,75.15,79.06,80.34,88.3,30.3,62.09,44.24,34,42.28,37.18,97.91,87.61,46.95,65.33,98.98,59.82,13.1,79.26,38.2,96.46,98.93,11.7,30.91,35.49,13.26,86.04,29.63,14.88,4.67,23.5,33.85,33.48,64.52,78.42,5.58,33.68,49.43,98.86,76.62,68.75,7.06,10.32)
stop_position_sorted=sort(stop_position)


#Graphs
#Random segmentations and calculation of the effective size(100 000 simulations)

Nes=c()
for(j in 1:100000){
  random_segmentation=runif(length(stop_position),0,100)           #random segmentation of our genes X times (X = number of observed stop codons)
  random_segmentation_sorted=sort(random_segmentation)            
  random_segmentation_sorted=append(0,random_segmentation_sorted)
  random_segmentation_sorted=append(random_segmentation_sorted,100)
  length_random_segments=c()
  i=1
  while(i < length(random_segmentation_sorted)){
    length_random_segments = append(length_random_segments, random_segmentation_sorted[i+1]-random_segmentation_sorted[i])
    i=i+1
  }
  Relative_length=length_random_segments/sum(length_random_segments)    #Calculation of the simulated effective segments sizes
  Ne=1/sum(Relative_length^2)
  Nes=append(Nes, Ne)                                                  #Regroupment of simulations results
}


# Calculation of the effective segmentation sizes created by observed STOPs codons

stop_position_sorted=append(0, stop_position_sorted)
stop_position_sorted=append(stop_position_sorted,100)

length_real_segments=c()
i=1
while(i < length(stop_position_sorted)){
  length_real_segments = append(length_real_segments, stop_position_sorted[i+1]-stop_position_sorted[i])  
  i=i+1
}

relative_length=length_real_segments/sum(length_real_segments)
ne=1/sum(relative_length^2)

line_y=c(1:19000)      #Creation of line vector for graphical representation
line_x=rep(ne, 19000)


##Plot stopcodons on a segment from 0 to 100%


plot(x=stop_position, y=rep(1, length(stop_position)), ylim=c(0.8,1.2), pch=16, col="red", axes=FALSE, ylab="", xlab="")
lines(x=0:100, y=rep(1,101))
lines(x=c(0,0), y=c(0.98, 1.02))
lines(x=c(100,100), y=c(0.98, 1.02))
text(x=99.5, y=0.97, "100%")
text(x=0.5, y=0.97, "0%")


#Creation of the graph S1

pdf("SupplementalFig1.pdf", width=11, height=8)
par(fig=c(0,1,0,0.7))
hist(Nes, ylim = c(0,25000), main="Effective segmentation size observed with random segmentations (100 000 simulations)", xlab = "Effective segmentation size", ylab = "Number of observation")
lines(line_x, line_y, col="red")
legend("topleft", "Observed segementation effective size (Stop codons)", col="red", lty=c(1,1), box.lwd=0)
par(fig=c(0,1,0.7,1), new=T)
plot(x=stop_position, y=rep(1, length(stop_position)), ylim=c(0.8,1.2), pch=16, col="red", axes=FALSE, ylab="", xlab="", main="Position of Stop codon on CDS length (%)")
lines(x=0:100, y=rep(1,101))
lines(x=c(0,0), y=c(0.98, 1.02))
lines(x=c(100,100), y=c(0.98, 1.02))
text(x=99.5, y=0.92, "100%", cex=1)
text(x=0.5, y=0.92, "0%", cex=1)
dev.off()

############################################################## Frameshifts ############################################################## 

##Exact same code but stop codons positions are replaces by frameshifts


#Frameshifts positions
FS_position=c(48.39,80.95,84.32,27.84,26.94,65.35,56.58,13.36,17.98,32.11,33.53,35.31,44.39,15.05,84.88,4.86,75.69,97.23,1.8,13.22,14.99,54.78,53.73,13,66.26,95,27.43,90.95,31.17,20.08,88.42,95.27,84.4,59.65,40.64,33.63,18.05,14.25,35.8,98.22,3.1,10.89,49.57,48.18,46.77,10.46,26.67,64.07,70.5,48.31,51.4,15.79,49.39,28.41,29.58,82.87,25.97,21.85,26.7,53.7,73.49,41.385,8.24,4.63,1.22,33.10,69.27,32.87,57.56,89.96)


frameshift_sorted=sort(FS_position)

#Random segmentations and calculatio of the effective size

Nes=c()
for(j in 1:100000){
  random_segmentation=runif(length(frameshift_sorted),0,100)
  random_segmentation_sorted=sort(random_segmentation)
  random_segmentation_sorted=append(0,random_segmentation_sorted)
  random_segmentation_sorted=append(random_segmentation_sorted,100)
  length_random_segments=c()
  j=1
  while(j < length(random_segmentation_sorted)){
    length_random_segments = append(length_random_segments, random_segmentation_sorted[j+1]-random_segmentation_sorted[j])
    j=j+1
  }
  Relative_length=length_random_segments/sum(length_random_segments)
  Ne=1/sum(Relative_length^2)
  Nes=append(Nes, Ne)
}


#Calculation of the effective sizes observed by the segmentation with our frameshifts

frameshift_sorted=append(0, frameshift_sorted)
frameshift_sorted=append(frameshift_sorted,100)

length_real_segments=c()
i=1
while(i < length(frameshift_sorted)){
  length_real_segments = append(length_real_segments, frameshift_sorted[i+1]-frameshift_sorted[i])
  i=i+1
}

relative_length=length_real_segments/sum(length_real_segments)
ne=1/sum(relative_length^2)

line_y=c(1:20000)
line_x=rep(ne, 20000)


# Graphic S2

pdf("SupplementalFig2.pdf", width=11, height=8)
par(fig=c(0,1,0,0.7))
hist(Nes, ylim = c(0,25000), main="Effective segmentation size observed with random segmentations (100 000 simulations)", xlab = "Effective segmentation size", ylab = "Number of observation")
lines(line_x, line_y, col="red")
legend("topleft", "Observed segementation effective size (Frameshifts)", col="red", lty=c(1,1), box.lwd=0)
par(fig=c(0,1,0.7,1), new=T)
plot(x=FS_position, y=rep(1, length(FS_position)), ylim=c(0.8,1.05), pch=16, col="blue", axes=FALSE, ylab="", xlab="", main="Position of frameshifts on CDS length (%)")
lines(x=0:100, y=rep(1,101))
lines(x=c(0,0), y=c(0.98, 1.02))
lines(x=c(100,100), y=c(0.98, 1.02))
text(x=99.5, y=0.92, "100%", cex=1)
text(x=0.5, y=0.92, "0%", cex=1)
dev.off()

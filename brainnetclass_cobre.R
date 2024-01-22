rm(list=ls())
#install.packages("devtools")
library(devtools)
#install_github("jesusdaniel/graphclass")
library(graphclass)
#install.packages("R.matlab")
library(R.matlab)
#install.packages("remotes")
#remotes::install_github("sidchop/brainconn")
#remotes::install_github("sidchop/brainconn")
#devtools::install_github("sidchop/brainconn")
#install.packages("rlang")  # Install the latest version of rlang
library(ggplot2)
library(brainconn)
#vignette("brainconn")
#install.packages("ggseg3d")
library(ggseg3d)
library(igraph)
#if (!require("NetworkToolbox", character.only=T, quietly=T)) {
#  devtools::insta3ll_github("AlexChristensen/NetworkToolbox")
#}
library(NetworkToolbox)
################################################################################
data(COBRE.data)

Aa=COBRE.data$subject.label
Bb=COBRE.data$Y.cobre
Cc=COBRE.data$X.cobre
COBRE_FC <- list(
  SubID = Aa,
  Lable = Bb,
  COBREFC = Cc)
# Save the data structure in a .mat file
writeMat(con="COBRE_FC.mat", myTable = COBRE_FC)

# do a loop with chatgpt
subjID <- Aa
Labels <- Bb  # 1 refers to SZ 
FC <- Cc

FC_tri_hc <- array(0, dim = c(263, 263, 124))  # HC 70
FC_tri_sz <- array(0, dim = c(263, 263, 124))  # SZ 54

for (i in 1:length(subjID)) {
  cat("Processing subject", i, "\n")
  if (Labels[i] == 1) {
      FC_vec_sz <- FC[i, , drop = FALSE]
      FC_tri_sz[, , i] <- get_matrix(FC_vec_sz)
  } else {
      FC_vec_hc <- FC[i, , drop = FALSE]
      FC_tri_hc[, , i] <- get_matrix(FC_vec_hc)
  }
}

non_zero_indices <- apply(FC_tri_sz, 3, function(x) any(x != 0))
FC_tri_sz_filtered <- FC_tri_sz[, , non_zero_indices]
non_zero_indices <- apply(FC_tri_hc, 3, function(x) any(x != 0))
FC_tri_hc_filtered <- FC_tri_hc[, , non_zero_indices]
#writeMat(con="COBRE_FC_FIN_Subjs.mat", HC=FC_tri_sz_filtered, SZ=FC_tri_hc_filtered)


FC_HC_avg <- apply(FC_tri_hc_filtered, c(1, 2), mean)  # dim(FC_HC_avg)
FC_SZ_avg <- apply(FC_tri_sz_filtered, c(1, 2), mean)

writeMat(con="COBRE_FC_FIN.mat", HC=FC_HC_avg, SZ=FC_SZ_avg)

data(power.parcellation)
# Node assignments (note that node 75 is missing on COBRE)
node.assignments <- power.parcellation$Master.Assignment[-75]
communities = lapply(c(1:13, -1), function(x) which(node.assignments==x))

Clabels=c("SMH","SMM","CO","AU","DM","MR","VS","FP","SA","SC","VA", "DA","CB","UN")
plot_adjmatrix(FC_HC_avg, communities = communities,community_labels = Clabels, colorlims =c(-0.2, 0.6), axislabel="brain regions", 
               main='HC')

plot_adjmatrix(FC_SZ_avg, communities = communities,community_labels = Clabels, colorlims =c(-0.2, 0.6), axislabel="brain regions", 
               main='SZ')

load("/home/liqiang/CVP/Neuroimaging/ISBI2024/power264.rda")
threshold_value <- 0.65
FC_HC_avg_thresholded <- ifelse(FC_HC_avg < threshold_value, 0, FC_HC_avg)
threshold_value <- 0.65
FC_SZ_avg_thresholded <- ifelse(FC_SZ_avg < threshold_value, 0, FC_SZ_avg)

#list_atlases()
brainconn(atlas ="power264", conmat=FC_HC_avg_thresholded, node.size=6, edge.alpha = 0.8, 
          edge.color.weighted = TRUE,   edge.width = 2, view="ortho",
          labels = FALSE)
#brainconn3D(atlas ="power264", conmat=FC_HC_avg_thresholded, show.legend = F)
brainconn(atlas ="power264", conmat=FC_SZ_avg_thresholded, node.size=6, edge.alpha = 0.8, 
          edge.color.weighted = TRUE,   edge.width = 2, view="ortho",
          labels = FALSE)
#Difference
FC_DC=FC_HC_avg-FC_SZ_avg
plot_adjmatrix(FC_DC, communities = communities, community_labels = Clabels,  colorlims =c(-0.2, 0.2), 
               axislabel="brain regions",  main='HC-SZ')

threshold_value <- 0.16
FC_DC_thresholded_strong <- ifelse(FC_DC < threshold_value, 0, FC_DC)

threshold_value <- -0.16
FC_DC_thresholded_weak <- ifelse(FC_DC > threshold_value, 0, FC_DC)

brainconn(atlas ="power264", conmat=FC_DC_thresholded_strong, node.size=6, edge.alpha = 0.8, 
          edge.color.weighted = TRUE,   edge.width = 2, view="ortho",
          labels = FALSE)

brainconn(atlas ="power264", conmat=FC_DC_thresholded_weak, node.size=6, edge.alpha = 0.8, 
          edge.color.weighted = TRUE,   edge.width = 2, view="ortho",
          labels = FALSE)

########################### graph
# read parcel labels for each subject
parcel.comm.path <- "/home/liqiang/CVP/Neuroimaging/ISBI2024/parcel_community"
parlabel <- data.frame(parcel_num=34453, 
                       community=NA, 
                       comm_label=NA, 
                       comm_shortlabel=NA)
plotlabel <- read.csv("/home/liqiang/CVP/Neuroimaging/ISBI2024/systemlabel_cobrecmi.txt", header=F,
                      col.names = c("community","comm_label","color","comm_shortlabel"))

# DESCRIPTION:
#   The degree to which edges are more dense within communites and more 
#   sparse between communities, which quantifies the segregation of a 
#   weighted network (Chan et al. 2014). 
#  
#  Inputs:   M,		Correlation matrix 
#            Ci,	Community affiliation vector (e.g., system labels)
#  Optional: 
#		        diagzero, 	Booleen for setting diagonal of input matrix to 0. 
#					              Default=TRUE
#		        negzero, 	  Booleen for setting negative edges of input matrix 
#					              to 0. Default=TRUE
#
#  Outputs: 
#     A list object is returned and contain the following elements - 
#           S,    System segregation calcualted with W & B                 
#           W,    Mean correlation between nodes within the same community     
#           B,    Mean correlation between nodes from different community   
# #########################################################################
#   Reference: Chan et al. (2014) PNAS E4997
#   Micaela Chan, UTD
#
#   Modification History:
#   May 2018: original (MYC)
#   Oct 2018: commented script (MYC)
########################################################################### 
segregation <- function(M=NULL, Ci=NULL, diagzero=TRUE, negzero=TRUE) {
  
  if(!isSymmetric.matrix(M, check.attributes = FALSE)) stop('Input matrix must be symmetric.')
  if(dim(M)[1]!=length(Ci)) stop('Length of community vector does not match matrix dimension.')
  
  # Set diagonal to 0
  if(diagzero==TRUE){
    diag(M) <- 0
  }
  
  # Set negatives to 0
  if(negzero==TRUE){
    M[which(M<0)] <- 0
  }
  
  within <- vector(mode = "numeric")
  between <- vector(mode = "numeric")
  
  Ci_order <- sort(unique(Ci))
  segresult <- list()
  
  for(i in 1:length(Ci_order)){ # loop through communities
    network <- Ci_order[i]
    
    ww = M[Ci==network,Ci==network]		# find index for within communitiy edges
    bb = M[Ci==network,which(Ci!=network)]	# find index for within communitiy edges
    
    within <- append(within, ww[upper.tri(ww)])
    between <- append(between, as.vector(bb))  
  }
  
  segresult$W <- mean(within)	# mean within community edges
  segresult$B <- mean(between)	# mean between community edges
  segresult$S <- (segresult$W-segresult$B)/segresult$W # calculate system segregation
  
  return(segresult)
}

Power14indsr <- c(
  13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 254,
  42, 43, 44, 45, 46,
  47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
  61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,
  74, 75, 76, 77, 78, 79, 80, 81, 82, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 136, 138,
  132, 133, 134, 135, 220,
  142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172,
  173, 174, 175, 176, 177, 178, 179, 180, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201,
  202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219,
  221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233,
  137, 234, 235, 236, 237, 238, 239, 240, 241,
  250, 251, 255, 256, 257, 258, 259, 260, 261, 262, 263,
  242, 243, 244, 245,
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 83, 84, 131, 139, 140, 141, 181, 182, 183, 184, 246, 247, 248, 249, 252, 253
)

Power14lbls <- c(
  rep(1, 30), rep(2, 5), rep(3, 14), rep(4, 13), rep(5, 57), rep(6, 5), rep(7, 31), rep(8, 25), rep(9, 18), rep(10, 13),
  rep(11, 9), rep(12, 11), rep(13, 4), rep(14, 28)
)

Power14indsr <- matrix(Power14indsr, nrow = 1)
Power14indsr <- t(Power14indsr)

#segH=segregation(FC_HC_avg, Ci=Power14lbls, diagzero=TRUE, negzero=TRUE)
#segS=segregation(FC_SZ_avg, Ci=Power14lbls, diagzero=TRUE, negzero=TRUE)

SegH_W<- numeric(70)
SegH_B<- numeric(70)
SegH_S<- numeric(70)

SegS_W<- numeric(54)
SegS_B<- numeric(54)
SegS_S<- numeric(54)

for (i in 1:70) {
  cat("Processing subject", i, "\n")
  segHL=segregation(FC_tri_hc_filtered[, , i], Ci=Power14lbls, diagzero=TRUE, negzero=TRUE)
  SegH_W[i] <- segHL$W
  SegH_B[i] <- segHL$B
  SegH_S[i] <- segHL$S
}

for (i in 1:54) {
  cat("Processing subject", i, "\n")
  segSL=segregation(FC_tri_sz_filtered[, , i], Ci=Power14lbls, diagzero=TRUE, negzero=TRUE)
  SegS_W[i] <- segSL$W
  SegS_B[i] <- segSL$B
  SegS_S[i] <- segSL$S
}

writeMat("COBRE_Seg.mat", SHW=SegH_W, SHB=SegH_B, SHS=SegH_S, SSW=SegS_W, SSB=SegS_B, SSS=SegS_S)

################################################################################
data(UMich.data)

BbB=UMich.data$Y.umich
CcC=UMich.data$X.umich

# do a loop with chatgpt
Labels <- BbB  # 1 refers to SZ 
FC <- CcC

FC_tri_hc <- array(0, dim = c(263, 263, 40))  # HC 40
FC_tri_sz <- array(0, dim = c(263, 263, 39))  # SZ 39

for (i in 1:79) {
  cat("Processing subject", i, "\n")
  if (Labels[i] == 1) {
    for (m in 1:40) {
      FC_vec_hc <- FC[i, , drop = FALSE]
      FC_tri_hc[, , m] <- get_matrix(FC_vec_hc, type = "undirected")
    }
  } else {
    for (n in 1:39) {
      FC_vec_sz <- FC[i, , drop = FALSE]
      FC_tri_sz[, , n] <- get_matrix(FC_vec_sz, type = "undirected")
    }
  }
}

FC_HC_avg <- apply(FC_tri_hc, c(1, 2), mean)  # dim(FC_HC_avg)
FC_SZ_avg <- apply(FC_tri_sz, c(1, 2), mean)
writeMat(con="Umich_FC_FIN.mat", HC=FC_HC_avg, SZ=FC_SZ_avg)

data(power.parcellation)
# Node assignments (note that node 75 is missing on COBRE)
node.assignments <- power.parcellation$Master.Assignment[-75]
communities = lapply(c(1:13, -1), function(x) which(node.assignments==x))

Clabels=c("SMH","SMM","CO","AU","DM","MR","VS","FP","SA","SC","VA", "DA","CB","UN")
plot_adjmatrix(FC_HC_avg, communities = communities,community_labels = Clabels, colorlims =c(-0.6, 0.8), axislabel="brain regions", 
               main='HC')

plot_adjmatrix(FC_SZ_avg, communities = communities,community_labels = Clabels, colorlims =c(-0.6, 0.8), axislabel="brain regions", 
               main='SZ')

load("/home/liqiang/CVP/Neuroimaging/ISBI2024/power264.rda")
threshold_value <- 0.76
FC_HC_avg_thresholded <- ifelse(FC_HC_avg < threshold_value, 0, FC_HC_avg)
threshold_value <- 0.70
FC_SZ_avg_thresholded <- ifelse(FC_SZ_avg < threshold_value, 0, FC_SZ_avg)

#list_atlases()
brainconn(atlas ="power264", conmat=FC_HC_avg_thresholded, node.size=6, edge.alpha = 0.8, 
          edge.color.weighted = TRUE,   edge.width = 2, view="ortho",
          labels = FALSE)
#brainconn3D(atlas ="power264", conmat=FC_HC_avg_thresholded, show.legend = F)
brainconn(atlas ="power264", conmat=FC_SZ_avg_thresholded, node.size=6, edge.alpha = 0.8, 
          edge.color.weighted = TRUE,   edge.width = 2, view="ortho",
          labels = FALSE)
#Difference
FC_DC=FC_HC_avg-FC_SZ_avg
plot_adjmatrix(FC_DC, communities = communities, community_labels = Clabels,  colorlims =c(-0.6, 0.6), 
               axislabel="brain regions",  main='HC-SZ')

threshold_value <- 0.9
FC_DC_thresholded_strong <- ifelse(FC_DC < threshold_value, 0, FC_DC)

threshold_value <- -0.9
FC_DC_thresholded_weak <- ifelse(FC_DC < threshold_value, 0, FC_DC)

brainconn(atlas ="power264", conmat=FC_DC_thresholded_strong, node.size=6, edge.alpha = 0.8, 
          edge.color.weighted = TRUE,   edge.width = 2, view="ortho",
          labels = FALSE)

brainconn(atlas ="power264", conmat=FC_DC_thresholded_weak, node.size=6, edge.alpha = 0.8, 
          edge.color.weighted = TRUE,   edge.width = 2, view="ortho",
          labels = FALSE)

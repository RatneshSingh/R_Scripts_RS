
#install.packages("adegenet","parallel","ape")
## import genepop format data to geneind format.
setwd("/home/ratnesh.singh/fangfang_rnaseq/Diversity_analysis_in_Turf_population/Diversity_analysis_in_Turf_population/adegenet")
library(adegenet)

SNPdata<-read.PLINK("/home/ratnesh.singh/fangfang_rnaseq/Diversity_analysis_in_Turf_population/Diversity_analysis_in_Turf_population/adegenet/run_as_41_groups_recodeA.raw",
 map.file="/home/ratnesh.singh/fangfang_rnaseq/Diversity_analysis_in_Turf_population/Diversity_analysis_in_Turf_population/adegenet/run_as_41_groups.map",
 chunkSize=10000,
 parallel=T ## leaving it sets it true and gives error in mclapply 
 )

#rm(list=ls())
#list.files()

#Basic analysis
postscript("SNP_visualization_plot.ps")
glPlot(SNPdata,posi="topleft")
dev.off()


### plot distribution of allele frequency
alFreq<-glMean(SNPdata)
postscript("Second_Allele_Frequency_plot.ps")
hist(alFreq,proba=T,col="red",xlab="Allele frequencies",main="Distribution of second allele frequencies")
temp<-density(alFreq)
lines(temp$x,temp$y*1.8,lwd=3)
dev.off()

## plot biallelic frequencies
#alFreq<-glMean(SNPdata) ## calculated above
bialFreq<-c(alFreq,1-alFreq)
postscript("Bi_Allele_Frequency_plot.ps")
hist(bialFreq,proba=T,col="darkseagreen",xlab="Allele frequencies",
 main="Distribution of second allele frequencies",nclass=20)
temp<-density(bialFreq,bw=.05)
lines(temp$x,temp$y*2,lwd=3)
dev.off()


## map missing data across loci
head(  glNA  (SNPdata),  20)

temp<- density( glNA (SNPdata), bw = 10)
postscript("SNP_visualization_NAvaues_plot.ps")
plot(temp, type = "n" , xlab = "Position in the alignment" , 
     main = "Location of the missing values (NAs)" , 
     xlim = c ( 0 , 200000 ))

polygon( c (temp $ x, rev (temp $ x)), 
         c (temp $ y, rep ( 0 , length (temp $ x))), 
         col = transp ( "blue" , .3 ))
points( glNA (SNPdata), rep ( 0 , nLoc (SNPdata)), pch = "|" , col = "blue")
dev.off()
##computations of distances between individuals
##seploc is used to create a list of smaller objects (here, 10 blocks of 10,000 SNPs
x<-  seploc(SNPdata, n.block =10 , parallel = T)

##dist is used within a lapplyloop to compute pairwise distances between individuals for each block
lD<-lapply(x, function(e) dist(as.matrix(e)))
D<-Reduce("+", lD)  ##distance matrix is obtained by summing these


## draw NJ tree using distance matrix
library(ape)
#postscript("NJtree_SNPdata_simple.ps")
#plot( njs(D),type="fan")
#title ( "A simple NJ tree turf 41 ind SNP data")
#dev.off()

##Run pca
pca.res=glPca(SNPdata,nf=41,parallel=T)


## plot eigenvalues
postscript("PCA_eigenvalues_plot_SNP.ps")
barplot(pca.res$eig, main="eigenvalues", col=heat.colors(length(pca.res$eig)))
dev.off()

## draw scatter plot of result
postscript("PCA_scatter_plot_SNP.ps")
scatter(pca.res, posi="bottomright",type = c( "lines"))
title("PCA of turf grass SNP data from 41 varieties\n axes 1-2")
dev.off()

##simple neighbour-joining (NJ) tree using distance data :
library(ape)
tre<-nj(dist(as.matrix(SNPdata)))
tre
postscript("PCA_NJtree.ps")
plot(tre, typ = "fan" , cex = 0.7)
title(  "NJ tree of turf grass SNP data from 41 varieties")
dev.off()

##The correspondance between both analyses can be better assessed using colors based on PCs; this is achieved by colorplot
postscript("PCA_scatter_colorplot.ps")
myCol<-colorplot(pca.res$scores,pca.res$scores,transp=T, cex=4)
abline(h=0,v=0,col="grey")
add.scatter.eig(pca.res$eig[1:40],2, 1, 2, posi="bottomright" ,inset=.05, ratio=.2)
dev.off()

##color tree
postscript("PCA_NJ_tree_color.ps")
plot(tre, typ="fan", show.tip=T,cex = 0.7)
tiplabels(pch = 20, col =myCol, cex = 4)
title(  "NJ tree of turf grass SNP data from 41 varieties")
dev.off()


##Discriminant Analysis of Principal Components (DAPC)
dapc1<-dapc(SNPdata, n.pca=10 ,n.da=1,parallel=F)

##DAPC can still make some decent discrimination
postscript("DAPC_scatter_plot_SNP.ps")
scatter(dapc1, scree.da=FALSE ,bg="white" ,
        posi.pca="topright", legend=TRUE, 
        txt.leg= paste("group", 1:2),
        col=c("red","blue"))
dev.off()

##composition plot confrms that groups are not entirely disentangled.
postscript("DAP_Ccompplot_plot_SNP.ps")
compoplot(dapc1, col = c ( "red" , "blue" ), lab = "" , txt.leg = paste ( "group" , 1 : 2 ), ncol = 2 )
dev.off()

## loading plot identifes pretty well the most discriminating alleles
postscript("DAPC_loading_plot_SNP.ps")
loadingplot (dapc1 $ var.contr, thres = 1e-3 )
dev.off()

#zoom in to the contributions of the last 100 SNPs to make sure that the tail indeed corresponds to the 50 last structured loci
postscript("DAPC_loading_plot_Zoomwd.ps")
loadingplot ( tail (dapc1 $ var.contr[, 1 ], 100 ), thres = 1e-3 )
dev.off()

#save all the data as image
save.image("/home/ratnesh.singh/fangfang_rnaseq/Diversity_analysis_in_Turf_population/Diversity_analysis_in_Turf_population/adegenet/run_adegenet_onturf.Rdata")

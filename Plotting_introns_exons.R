setwd("/home/ratnesh.singh/Pineapple/PinappGenome_Illinois/Ac_CG_contents/introns")
#load("/home/ratnesh.singh/Pineapple/PinappGenome_Illinois/Ac_CG_contents/introns/Plotting_introns_exons.RData")
library(cluster)
library(fpc)
#library(TSclust)
#file="pineapple.v3.20140921.assembly.intron-exon-gc"



file_list=c(
  "Ptrichocarpa_210.intron-exon-gc", 
  "Athaliana_167.intron-exon-gc",     
  "Cpapaya_113.intron-exon-gc", 
  "Vvinifera_145.intron-exon-gc",
  "Nnucifera_MegaScaffold.intron-exon-gc",
  "Spirodela_polyrhiza_genome.intron-exon-gc",
  "Pequestis_1213.scafSeq.FG2_superscaffold.intron-exon-gc", 
  "PdactyKAsm30_r20101206.intron-exon-gc", 
  "Macuminata_genomic_v1.intron-exon-gc", 
  "pineapple.v3.20140921.assembly.intron-exon-gc",  
  "Sbicolor_v2.1_255.intron-exon-gc",  
  "Zmays_284_AGPv3.intron-exon-gc",
  "Sitalica_164_v2.intron-exon-gc",                   
  "Osativa_204_v7.0.intron-exon-gc",
  "Bdistachyon_283_assembly_v2.0.intron-exon-gc",
  "Atrichopoda_291_v1.0.intron-exon-gc", 
  "Smoellendorffii_91.intron-exon-gc",
  "Vcarteri_199.intron-exon-gc"            
)


### plot gene length and CDS length lengt

resln=600  ## image resolution for saving.
psize=0.8
tiff(paste(dfile,"Exon_intron_GC","tiff",sep="."),width = 1, height = 18, units = 'in', res = resln,compression = "lzw")
par(mfrow=c(18,1))
for(i in 1:length(file_list)){
  
  file<-file_list[i]
 # file<-file_list[11]  
  print(paste("Processing file",file,sep=" "))
  #########################
  perlensliding=paste(file,"perlensliding.table",sep=".")
  GCperexon=paste(file,"GCperexon.table",sep=".")
  GCperintron=paste(file,"GCperintron.table",sep=".")
  Lenperexon=paste(file,"Lenperexon.table",sep="")
  Lenperintron=paste(file,"Lenperintron.table",sep=".")
  perlen=paste(file,"perlen.table",sep=".")
  all=paste(file,"table",sep=".")
  
  #dfile=perlensliding
  #dfile=GCperexon
  #dfile=GCperintron
  #dfile=Lenperexon
  #dfile=Lenperintron
  #dfile=perlen
  dfile=all
  
  #### all file contains these columns
  ## [1] "NumIntrons"              "gene_length"             "CDS_length"              "Intron_length"          
  ## [5] "X.GC_exons"              "X.GC_introns"            "GC_each_exons"           "GC_each_introns"        
  ## [9] "len_each_exon"           "len_each_intron"         "GC_per_10percent_of_CDS"
  print(paste("Reading data from file",dfile,sep=" "))
  data=read.table(dfile,header=T,row.names=1)
  dat <- data 
  
  
  pars=colnames(data)
 
  lt60=data.frame()
  gt60=data.frame()
  ltlist=list()
  gtlist=list()
  for (i in 1:length(data[,1])){
    if(data[i,5] < 60){
      ltlist<-c(ltlist,i)
    }
    else{
      gtlist<-c(gtlist,i)
    }
  }
  
  gtlist<-as.numeric(gtlist)
  ltlist<-as.numeric(ltlist)
  
  gt60<-data[gtlist,]
  lt60<-data[ltlist,]
  ltcor<-cor(as.numeric(lt60[,5]),as.numeric(lt60[,6]))  
  gtcor<-cor(as.numeric(gt60[,5]),as.numeric(gt60[,6]))
  
  plot( lt60[,5],  lt60[,6],   col="blue",   xlab="% GC exons",   ylab="% GC introns",   ylim=c(0,100),   xlim=c(0,100),   sub=paste("Coeff(<60%):",  round(ltcor,2),"Coeff(>60%):",round(gtcor,2),sep="    "), main=file)
  points(gt60[,5],gt60[,6],col="red") 

  
 
}
dev.off()















##removed;            "Ppatens_152.intron-exon-gc", 
for(i in 1:length(file_list)){
  
  file<-file_list[i]
  print(paste("Processing file",file,sep=" "))
  #########################
  
  resln=600  ## image resolution for saving.
  psize=0.8
  
  
  #########################
  perlensliding=paste(file,"perlensliding.table",sep=".")
  GCperexon=paste(file,"GCperexon.table",sep=".")
  GCperintron=paste(file,"GCperintron.table",sep=".")
  Lenperexon=paste(file,"Lenperexon.table",sep="")
  Lenperintron=paste(file,"Lenperintron.table",sep=".")
  perlen=paste(file,"perlen.table",sep=".")
  all=paste(file,"table",sep=".")
  
  #dfile=perlensliding
  #dfile=GCperexon
  #dfile=GCperintron
  #dfile=Lenperexon
  #dfile=Lenperintron
  #dfile=perlen
  dfile=all
  
  #### all file contains these columns
  ## [1] "NumIntrons"              "gene_length"             "CDS_length"              "Intron_length"          
  ## [5] "X.GC_exons"              "X.GC_introns"            "GC_each_exons"           "GC_each_introns"        
  ## [9] "len_each_exon"           "len_each_intron"         "GC_per_10percent_of_CDS"
  print(paste("Reading data from file",dfile,sep=" "))
  data=read.table(dfile,header=T,row.names=1)
  dat <- data 
  
  
  pars=colnames(data)
  
  par1="NumIntrons"
  par2="gene_length"
  tiff(paste(dfile,"AllGraphsInOne","tiff",sep="."),width = 25, height = 25, units = 'in', res = resln,compression = "lzw")
  par(mfrow=c(5,5))
  for (i in 1:5){  
    for (j in 2:6){
      if (i == j) next 
      par1=pars[i]
      par2=pars[j]
      corcoef<-cor(as.numeric(data[,par1]),as.numeric(data[,par2]))
      
      ### save image as postscript file
      file=paste(dfile,par1,"vs",par2,"ps",sep=".")
      #print(file)
      #postscript(file)
      #plot(data[,par1],data[,par2],main=paste(par1,"vs",par2,sep=" "),xlab=par1,ylab=par2)
      #dev.off()
      ### save as tiff too
      file=paste(dfile,par1,"vs",par2,"tiff",sep=".")
      print(paste("Printing plot for:",file))
      #tiff(file,width = 7, height = 7, units = 'in', res = resln,compression = "lzw")
      plot(data[,par1],data[,par2],main=paste(par1,"vs",par2,sep=" "),xlab=par1,ylab=par2,sub=paste("Correlation coefficient:",round(corcoef,3)),cex=psize)
      #dev.off()
      
      
    }
  }
  dev.off()
}
  
  # Kmeans clustre analysis
  #clus <- kmeans(dat, centers=2,iter.max = 100)
  #plotcluster(dat, clus$cluster)
  
  #### process each exon/intron data and plot them.
  ## split each element with comma and create alist of unequal length 
  geev<-strsplit(as.character(data[,7]),",")
  
  ### add values to data.frame
  #geed<-data.frame()
  #for (i in seq(along=geev)){
  #  for (j in 1:as.numeric(length(geev[[i]]))){
  #       geed[i,j]<-geev[[i]][j]        
  
  #  }
  #}
  
  
  max.len <- max(sapply(geev, length))
  corrected.list <- lapply(geev, function(x) {c(x, rep(NA, max.len - length(x)))})
  geem <- do.call(rbind, corrected.list)

  ### plot the matrix rows
  tiff(paste(dfile,"GC_per_Exon","tiff",sep="."),width = 7, height = 7, units = 'in', res = resln,compression = "lzw")
  matplot(t(geem),type="p",pch=1,col=1:1,xlab="Exon Numbers",ylab="GC %",main="")
  dev.off()
  
  ### plot the last exons only
  first_exons=list()
  last_exons=list()
  num_exons=list()
  second_exons=list()
  sec_num_exons=list()
  third_exons=list()
  third_num_exons=list()
  for(i in 1:nrow(geem)){
    index=length(which(!is.na(geem[i,])))
    if(index == 1) {
      
      next
    }
    x<-geem[i,index]
    num_exons<-c(num_exons,index)
    last_exons<-c(last_exons,x)
    first_exons<-c(first_exons,geem[i,1])
    if(index > 2){
      second_exons<-c(second_exons,geem[i,2])
      sec_num_exons<-c(sec_num_exons,index)
      
    }
    if(index > 3){
      third_exons<-c(third_exons,geem[i,3])
      third_num_exons<-c(third_num_exons,index)
      
    }
  }
  tiff(paste(dfile,"GC_first_vs_last_Exon","tiff",sep="."),width = 7, height = 7, units = 'in', res = resln,compression = "lzw")
  plot(last_exons,first_exons, xlab="GC% (Last Exon)", ylab="GC% (First Exon)",xlim=c(0,100),ylim=c(0,100))
  dev.off()
  
  tiff(paste(dfile,"GC_first_vs_num_Exon","tiff",sep="."),width = 7, height = 7, units = 'in', res = resln,compression = "lzw")
  plot(first_exons,num_exons,xlab="GC% (First Exon)",ylab="Num exons",xlim=c(0,100))
  dev.off()
  
  tiff(paste(dfile,"GC_last_vs_num_Exon","tiff",sep="."),width = 7, height = 7, units = 'in', res = resln,compression = "lzw")
  plot(last_exons,num_exons,xlab="GC% (Last Exon)",ylab="Num exons",xlim=c(0,100))
  dev.off()
  
  tiff(paste(dfile,"GC_second_vs_num_Exon","tiff",sep="."),width = 7, height = 7, units = 'in', res = resln,compression = "lzw")
  plot(second_exons,sec_num_exons,xlab="GC% (Second Exon)",ylab="Num exons",xlim=c(0,100))
  dev.off()
  
  tiff(paste(dfile,"GC_third_vs_num_Exon","tiff",sep="."),width = 7, height = 7, units = 'in', res = resln,compression = "lzw")
  plot(third_exons,third_num_exons,xlab="GC% (Third Exon)",ylab="Num exons",xlim=c(0,100))
  dev.off()
  
  
  
  
  ### multi image in one
  psize=0.5
  tiff(paste(dfile,"GC_first_vs_last_Exon","tiff",sep="."),width = 7, height = 7, units = 'in', res = resln,compression = "lzw")
  
  par(mfrow=c(3,2))
  lfe<-cor(as.numeric(last_exons),as.numeric(first_exons))
  plot(last_exons,first_exons, xlab="GC% (Last Exon)", ylab="GC% (First Exon)",xlim=c(0,100),ylim=c(0,100),cex=psize,sub=paste("Correlation coefficient:",round(lfe,3)))
  
  fne<-cor(as.numeric(first_exons),as.numeric(num_exons))
  plot(first_exons,num_exons,xlab="GC% (First Exon)",ylab="Num exons",xlim=c(0,100),cex=psize,sub=paste("Correlation coefficient:",round(fne,3)))
  
  lne<-cor(as.numeric(last_exons),as.numeric(num_exons))
  plot(last_exons,num_exons,xlab="GC% (Last Exon)",ylab="Num exons",xlim=c(0,100),cex=psize,sub=paste("Correlation coefficient:",round(lne,3)))
  
  sne<-cor(as.numeric(second_exons),as.numeric(sec_num_exons))
  plot(second_exons,sec_num_exons,xlab="GC% (Second Exon)",ylab="Num exons",xlim=c(0,100),cex=psize,sub=paste("Correlation coefficient:",round(sne,3)))
  
  
  tne<-cor(as.numeric(third_exons),as.numeric(third_num_exons))
  plot(third_exons,third_num_exons,xlab="GC% (Third Exon)",ylab="Num exons",xlim=c(0,100),cex=psize,sub=paste("Correlation coefficient:",round(tne,3)))
  
  
  dev.off()
  
  
  
  
  
  save.image("/home/ratnesh.singh/Pineapple/PinappGenome_Illinois/Ac_CG_contents/introns/Plotting_introns_exons.RData")
}

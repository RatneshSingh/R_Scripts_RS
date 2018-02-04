#setwd("/home/ratnesh.singh/Pineapple/PinappGenome_Illinois/Ac_CG_contents/introns")
setwd("C:\\Users\\ratnesh.singh\\Google Drive\\Ac_GC_Content analysis\\introns")
#load("/home/ratnesh.singh/Pineapple/PinappGenome_Illinois/Ac_CG_contents/introns/Plotting_introns_exons.RData")
library(cluster)
library(fpc)

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

name_list=c(
  "P.trichocarpa", 
  "A.thaliana",     
  "C.papaya", 
  "V.vinifera",
  "N.nucifera",
  "Spirodela polyrhiza",
  "P.equestis", 
  "P.dactylifera", 
  "M.acuminata", 
  "pineapple",  
  "S.bicolor",  
  "Z.mays",
  "S.italica",                   
  "O.sativa",
  "B.distachyon",
  "A.trichopoda", 
  "S.moellendorffii",
  "V.carteri"            
)





### plot CDS GC and Intron GC contents.

resln=600  ## image resolution for saving.
psize=0.8
tiff(paste("GC_exon_lengths_intron_all_samples","tiff",sep="."),width = 20, height = 25, units = 'in', res = resln,compression = "lzw")
par(mfrow=c(5,4))
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
  ## [1] "NumIntrons"              [2]"gene_length"             [3]"CDS_length"              [4]"Intron_length"          
  ## [5] "X.GC_exons"              [6]"X.GC_introns"            [7]"GC_each_exons"           [8]"GC_each_introns"        
  ## [9] "len_each_exon"           [10]"len_each_intron"         [10]"GC_per_10percent_of_CDS"
  
  print(paste("Reading data from file",dfile,sep=" "))
  data=read.table(dfile,header=T,row.names=1)
  pars=colnames(data)
  pars[5]<-"% GC exons"
  pars[6]<-"% GC introns"
  ## compare columns
 xc=3
 yc=5
 
 
 ### filter data based on GC exon ###
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
############################################################  
  gt60<-data[gtlist,]
  lt60<-data[ltlist,]
  ltcor<-cor(as.numeric(lt60[,xc]),as.numeric(lt60[,yc]))  
  gtcor<-cor(as.numeric(gt60[,xc]),as.numeric(gt60[,yc]))

  maxlen<-max(c(as.numeric(lt60[,xc]),as.numeric(lt60[,yc])))
 
#### Asign plotting data to correct variables  ##### 
 x1<-lt60[,xc]
 y1<-lt60[,yc]

 x2<-gt60[,xc]
 y2<-gt60[,yc]


 xlabl<-"CDS length"
 ylabl<-"% GC introns"

 xlimt<-c(0,maxlen)
 ylimt<-c(0,100)
 subt<-paste("Coeff(<60%):",  round(ltcor,2),"Coeff(>60%):",round(gtcor,2),sep="    ")
 maint<-name_list[i]
 
 #### plot the data ####
  
  plot( x1,  y1,   col="blue",   xlab=xlabl,   ylab=ylabl,   ylim=ylimt,   xlim=xlimt,   sub=subt, main=maint)
  points(x2, y2,   col="red") 

  
 
}
dev.off()

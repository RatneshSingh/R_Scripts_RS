#### all file contains these columns. provide which two columns you want to compare
## [1] "NumIntrons"              [2]"gene_length"             [3]"CDS_length"              [4]"Intron_length"          
## [5] "X.GC_exons"              [6]"X.GC_introns"            [7]"GC_each_exons"           [8]"GC_each_introns"        
## [9] "len_each_exon"           [10]"len_each_intron"         [11]"GC_per_10percent_of_CDS"
pars<-c("NumIntrons","gene_length","CDS_length","Intron_length","X.GC_exons","X.GC_introns",
        "GC_each_exons","GC_each_introns","len_each_exon" ,"len_each_intron","GC_per_10percent_of_CDS")
xc=5  ## plot this column on x axis
yc=6  ## plot this column on y axis
div_col=4  ## divide the data in to parts based on column.
div_val=1000  ## divide the data in to parts based on cutoff.








########################################################################################################
#### Dont change anything below. if files list is not correct, changes can be done.
########################################################################################################
print(paste("Calculating and plotting",pars[xc],pars[yc],sep=" "))
setwd("/home/ratnesh.singh/Pineapple/PinappGenome_Illinois/Ac_CG_contents/introns")
#setwd("C:\\Users\\ratnesh.singh\\Google Drive\\Ac_GC_Content analysis\\introns")
#load("/home/ratnesh.singh/Pineapple/PinappGenome_Illinois/Ac_CG_contents/introns/Plotting_introns_exons.RData")


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
  "Pabies1.0-genome-gene-only.intron-exon-gc",
  "Smoellendorffii_91.intron-exon-gc",
  "Ppatens_251_v3.softmasked.intron-exon-gc",
  "Vcarteri_199.intron-exon-gc"            
 )

name_list=c(
  "Populus trichocarpa", 
  "Arabidopsis thaliana",     
  "Carica papaya", 
  "Vitis vinifera",
  "Nelumbo nucifera",
  "Spirodela polyrhiza",
  "Phalaenopsis equestis", 
  "Phoenix dactylifera", 
  "Musa acuminata", 
  "Ananas comosus",  
  "Sorghum bicolor",  
  "Zea mays",
  "Setaria italica",                   
  "Oryza sativa",
  "Brachypodium distachyon",
  "Amborella trichopoda",
  "Picea abies",
  "Selaginella moellendorffii",
  "Physcomitrella patens",
  "Volvox carteri"
)


#####################################################################################################
library(cluster)
library(fpc)

### plot CDS GC and Intron GC contents.

resln=600  ## image resolution for saving.
psize=0.8
outfile<-paste(pars[xc],pars[yc],length(file_list),"samples","colored_red",pars[div_col],"gt",div_val,"tiff",sep=".")
outfile=gsub(" ", "", outfile, fixed = TRUE)
tiff(outfile,width = 20, height = 25, units = 'in', res = resln,compression = "lzw")
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
  
  ##################################################################################################### 
  
  print(paste("Reading data from file",dfile,sep=" "))
  data=read.table(dfile,header=T,row.names=1)
  pars=colnames(data)
  pars[5]<-"% GC exons"
  pars[6]<-"% GC introns"
  pars[7]<-"% GC each exons"
  pars[8]<-"% GC each introns"
  pars[11]<-"% GC per 10percent of length"
  ## compare columns
  
  
  
  ### filter data based on GC exon ###
  lt60=data.frame()
  gt60=data.frame()
  ltlist=list()
  gtlist=list()
  for (j in 1:length(data[,1])){
    if(data[j,div_col] < div_val){
      ltlist<-c(ltlist,j)
    }
    else{
      gtlist<-c(gtlist,j)
    }
  }
  
  gtlist<-as.numeric(gtlist)
  ltlist<-as.numeric(ltlist)
  ############################################################  
  gt60<-data[gtlist,]
  lt60<-data[ltlist,]
  ltcor<-cor(as.numeric(lt60[,xc]),as.numeric(lt60[,yc]))  
  gtcor<-cor(as.numeric(gt60[,xc]),as.numeric(gt60[,yc]))
  
  xmax<-max(c(as.numeric(lt60[,xc]),as.numeric(gt60[,xc])))
  ymax<-max(c(as.numeric(lt60[,yc]),as.numeric(gt60[,yc])))
  
  #### Asign plotting data to correct variables  ##### 
  x1<-lt60[,xc]
  y1<-lt60[,yc]
  
  x2<-gt60[,xc]
  y2<-gt60[,yc]
  
  ### assign correct labels
  xlabl<-pars[xc]
  ylabl<-pars[yc]
  
  ### deciding what are axis values
  if(xc < 5 | xc > 8){
    xlimt<-c(0,xmax)
  }
  if(yc < 5 | yc > 8){
    ylimt<-c(0,ymax)
  }
  
  if(xc > 4 & xc < 9){
    xlimt<-c(0,100)
  }
  if(yc > 4 & yc < 9){
    ylimt<-c(0,100)
  }
  
  
  subt<-paste("Coeff(",pars[div_col],"<",div_val,"):",  round(ltcor,2),"    Coeff(",pars[div_col],">",div_val,"):",  round(gtcor,2),sep=" ")
  maint<- substitute(expr = italic(m),env = list(m=name_list[i]))
  
  #### plot the data ####
  
  plot( x1,  y1,   col="blue",   xlab=xlabl,   ylab=ylabl,   ylim=ylimt,   xlim=xlimt,   sub=subt, main=maint)
  points(x2, y2,   col="red") 
  
  
  
}
dev.off()

print(paste("Saved image as",outfile,sep=" "))

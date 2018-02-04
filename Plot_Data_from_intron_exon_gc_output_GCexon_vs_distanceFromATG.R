## args will be stored as arg[1] arg[2] etc..
#args<-commandArgs(TRUE)

#### fancy argument reader
library(getopt)
spec <- matrix(c(
        'pg1'     ,'p', 1, "integer", "-1: ignore this filter; 0 : plot less than div_val data only; 1 : plot greater than div_val data only; 2 : plot equal to div_val data only; 3: plot data from lt to gt; 4: plot all the data",
        'dc1'   ,'c', 1, "integer", "col numberto use for filter1.Values: [1]NumIntrons[2]gene_length[3]CDS_length[4]Intron_length[5]X.GC_exons[6]X.GC_introns",
        'odv1'  , 'u', 2, "integer", "lower Value to be used for filter1",
        'odv2'  , 'v', 2, "integer", "upper Value to be used for filter1",
        'pg2'     ,'q', 1, "integer", "-1: ignore this filter; 0 : plot less than div_val data only; 1 : plot greater than div_val data only; 2 : plot equal to div_val data only; 3: plot data from lt to gt; 4: plot all the data",
        'dc2'   ,'e', 1, "integer", "col numberto use for filter1.Values: [1]NumIntrons[2]gene_length[3]CDS_length[4]Intron_length[5]X.GC_exons[6]X.GC_introns",
        'tdv1'  , 'w', 2, "integer", "lower Value to be used for filter2",
        'tdv2'  , 'x', 2, "integer", "upper Value to be used for filter2",
        'nintr' , 'i', 2, "integer", "Number of intron to plot. default:26",
        'coeff' , 'f', 2, "numerical", "Coefficient cutoff to color red(positive correlation)/blue(negative correlation)",
        'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);


if (!is.null(opt$help) || is.null(opt$pg1)) {
    cat(paste(getopt(spec, usage=T),"\n"));
    q();
}

## default values if not given.
if(is.null(opt$pg1)) opt$pg1=4
if(is.null(opt$dc1))opt$dc1=1
if(is.null(opt$odv1)) opt$odv1=26
#if(is.null(opt$odv2)) opt$odv2=26

if(is.null(opt$pg2)) opt$pg2=-1
if(is.null(opt$dc2)) opt$dc2=5
#if(is.null(opt$tdv1)) opt$tdv1=0
#if(is.null(opt$tdv2)) opt$tdv2=60
### create div value for filter
if(is.null(opt$nintr))opt$nintr=26
if(is.null(opt$coeff)) opt$coeff=0.5
#### all the arguments will be in the form of opt$argname

#### all file contains these columns. provide which two columns you want to compare
## [1] "NumIntrons"              [2]"gene_length"             [3]"CDS_length"              [4]"Intron_length"          
## [5] "X.GC_exons"              [6]"X.GC_introns"            [7]"GC_each_exons"           [8]"GC_each_introns"        
## [9] "len_each_exon"           [10]"len_each_intron"         [11]"GC_per_10percent_of_CDS"
pars<-c("NumIntrons","gene_length","CDS_length","Intron_length","X.GC_exons","X.GC_introns",
        "GC_each_exons","GC_each_introns","len_each_exon" ,"len_each_intron","GC_per_10percent_of_CDS")
#xc=5  ## plot this column on x axis
#yc=4  ## plot this column on y axis
div_col=opt$dc1  ## divide the data in to parts based on column.

if(!is.null(opt$odv1)) {div_val=c(opt$odv1)
}else if(!is.null(opt$odv1) && !is.null(opt$odv2)) {div_val=c(opt$odv1,opt$odv2)}   ## cutoff value for dividing the data
plot_gt=opt$pg1  ## 0 : plot less than div_val data only; 1 : plot greater than div_val data only; 2 : plot equal to div_val data only; 3: plot data from lt to gt; 4: plot all the data.

if(plot_gt==0) title_gt="lt"
if(plot_gt==1) title_gt="gt"
if(plot_gt==2) title_gt="eq"
if(plot_gt==3) title_gt="lt2gt" 
if(plot_gt==4) title_gt="ltNgt"

if(plot_gt == 3){if(length(div_val) < 2 ){plot_gt=4}}
if(plot_gt == 4){div_col=0;div_val=0}

########################################################################################################################################################
###Second filter
div_col2=opt$dc2  ## divide the data in to parts based on column.
if(!is.null(opt$tdv1)) {div_val2=c(opt$tdv1)
}else if(!is.null(opt$tdv1) && !is.null(opt$tdv2)) {div_val2=c(opt$tdv1,opt$tdv2)}  ## cutoff value for dividing the data
plot_gt2=opt$pg2  ## -1: No filtering; 0 : plot less than div_val data only; 1 : plot greater than div_val data only; 2 : plot equal to div_val data only; 3 : plot all the data.

if(plot_gt2 < 0) title_gt2=""
if(plot_gt2==0) title_gt2="lt"
if(plot_gt2==1) title_gt2="gt"
if(plot_gt2==2) title_gt2="eq"
if(plot_gt2==3) title_gt2="lt2gt"
if(plot_gt2==4) title_gt2="All"
if(plot_gt2 == 3){if(length(div_val2) < 2 ){plot_gt2=4}}
if(plot_gt2 == 4){div_col2=0;div_val2=0}


### what to plot on which axis 
##'gee':gc_each-exon
##'gei':gc_each_intron
##'exd':dist_each-exon
##'ind':dist_each_intron
on_x = 'gee'
on_y = 'gei'

if(on_x == 'gee') labx="Exon number"
if(on_x == 'gei') labx="Intron number"
if(on_x == 'exd') labx="Exon number"
if(on_x == 'ind') labx="Intron number"
if(on_y == 'gee') laby="Exon number"
if(on_y == 'gei') laby="Intron number"
if(on_y == 'exd') laby="Exon number"
if(on_y == 'ind') laby="Intron number"




Num_introns=opt$nintr  ## Number of introns to plot for.
coeff_cutoff=opt$coeff   ### correlation cutoff for coloring red




title1=""
title2=""
if(plot_gt  > 0) {title1<-paste("filt1On",pars[div_col],title_gt,paste(div_val,collapse="_"),sep=".")}
if(plot_gt2 > 0) {title2<-paste("fil2On", pars[div_col2],title_gt2,paste(div_val2,collapse="_"),sep=".")}

########################################################################################################
#### Dont change anything below. if files list is not correct, changes can be done.
########################################################################################################
#print(paste("Calculating and plotting",pars[xc],pars[yc],sep=" "))
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
#coeff_cutoff=0.3  ### for colring the red in plot.
outfile<-paste(on_x,"vs",on_y,title1,title2,"NumIntronsShown",Num_introns,length(file_list),"samples","colored_red.Corr_coeff","gt",coeff_cutoff,"tiff",sep=".")
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
  
  ####### collect information for columns 7,8,9,10,11 and tranlate to matrix by seperating to comma.
if(on_x == 'gee' || on_y == 'gee'){  
  gee<-strsplit(as.character(data[,7]),",")
  gee<-lapply(gee,function(x){as.numeric(x)})
  max.len <- max(sapply(gee, length))
  corrected.list <- lapply(gee, function(x) {c(as.numeric(x), rep(NA, max.len - length(x)))})
  geem <- do.call(rbind, corrected.list)
}
if(on_x == 'gei' || on_y == 'gei'){  
  gei<-strsplit(as.character(data[,8]),",")
  gei<-lapply(gei,function(x){as.numeric(x)})
  max.len <- max(sapply(gei, length))
  corrected.list <- lapply(gei, function(x) {c(as.numeric(x), rep(NA, max.len - length(x)))})
  geim <- do.call(rbind, corrected.list)
}
if(on_x == 'lee' || on_y == 'lee'){  
  lee<-strsplit(as.character(data[,9]),",")
  lee<-lapply(lee,function(x){as.numeric(x)})
  max.len <- max(sapply(lee, length))
  corrected.list <- lapply(lee, function(x) {c(as.numeric(x), rep(NA, max.len - length(x)))})
  leem <- do.call(rbind, corrected.list)
}
if(on_x == 'lei' || on_y == 'lei'){  
  lei<-strsplit(as.character(data[,10]),",")
  lei<-lapply(lei,function(x){as.numeric(x)})
  max.len <- max(sapply(lei, length))
  corrected.list <- lapply(lei, function(x) {c(as.numeric(x), rep(NA, max.len - length(x)))})
  leim <- do.call(rbind, corrected.list)
}
if(on_x == 'gpdec' || on_y == 'gpdec'){  
  gpdec<-strsplit(as.character(data[,11]),",")
  gpdec<-lapply(gpdec,function(x){as.numeric(x)})
  max.len <- max(sapply(gpdec, length))
  corrected.list <- lapply(gpdec, function(x) {c(as.numeric(x), rep(NA, max.len - length(x)))})
  gpdecm <- do.call(rbind, corrected.list)
}  
  ##### calculate distances of each exon and intron from ATG. Convert them into matrix
if(on_x == 'exd' || on_y == 'exd'){  
  dist_exo <- lapply(seq_along(lee),function(ei) {lapply(1:length(lee[[ei]]),function(ej){sum(lee[[ei]][0:(ej-1)],lei[[ei]][0:(ej-1)])})})
  max.len <- max(sapply(dist_exo, length))
  corrected.list <- lapply(dist_exo, function(x) {c(as.numeric(x), rep(NA, max.len - length(x)))})
  dist_exom <- do.call(rbind, corrected.list)
}
if(on_x == 'ind' || on_y == 'ind'){  
  dist_int <- lapply(seq_along(lei),function(ii) {lapply(1:length(lei[[ii]]),function(ij){sum(lee[[ii]][0:ij],lei[[ii]][0:(ij-1)])} )})
  max.len <- max(sapply(dist_int, length))
  corrected.list <- lapply(dist_int, function(x) {c(as.numeric(x), rep(NA, max.len - length(x)))})
  dist_intm <- do.call(rbind, corrected.list)
}
  ####################################################################################################
  ### group filter data based on other criteria ###
  ltlist=list()
  gtlist=list()
  eqlist=list()
  lt2gtlist=list()
  if(plot_gt < 4){
    for (j in 1:length(data[,1])){
      if     (plot_gt == 0 && data[j,div_col]  < div_val[[1]] ){ltlist<-c(ltlist,j) 
      }else if(plot_gt == 1 && data[j,div_col]  > div_val[[1]] ){gtlist<-c(gtlist,j)
      }else if(plot_gt == 2 && data[j,div_col] == div_val[[1]] ){eqlist<-c(eqlist,j) 
      }else if(plot_gt == 3 && length(div_val) > 1 && data[j,div_col] >= div_val[[1]] && data[j,div_col] <= div_val[[2]] ){lt2gtlist<-c(lt2gtlist,j) }
    }
  } 
     if(plot_gt == 0) {fil1_list<-ltlist
     }else if(plot_gt == 1) {fil1_list<-gtlist  
     }else if(plot_gt == 2) {fil1_list<-eqlist  
     }else if(plot_gt == 3){fil1_list<-lt2gtlist
     }else                  {fil1_list<-1:length(data[,1])}

# #### trying to create a second filter for multi level filtering
## new data set after above filtering.
ltlist=list()
gtlist=list()
eqlist=list()
lt2gtlist=list()
if(plot_gt2 < 4 && plot_gt2 > 0){
  for (j2 in fil1_list){
    if      (plot_gt2 == 0 && data[j2,div_col2]  < div_val2[[1]] ){ltlist<-c(ltlist,j2) 
    }else if(plot_gt2 == 1 && data[j2,div_col2]  > div_val2[[1]] ){gtlist<-c(gtlist,j2) 
    }else if(plot_gt2 == 2 && data[j2,div_col2] == div_val2[[1]] ){eqlist<-c(eqlist,j2) 
    }else if(plot_gt2 == 3 && length(div_val2) > 1 && data[j2,div_col2] >= div_val2[[1]] && data[j2,div_col2] <= div_val2[[2]] ){lt2gtlist<-c(lt2gtlist,j2) }
  }
} 

     if(plot_gt2 == 0){fil2_list<-ltlist
     }else if(plot_gt2 == 1){fil2_list<-gtlist
    }else if(plot_gt2 == 2){fil2_list<-eqlist
    }else if(plot_gt2 == 3){fil2_list<-lt2gtlist
    }else{fil2_list<- fil1_list}
###### end filtering multilevel and assign final filtered list to fil_list.
fil_list <-as.numeric(fil2_list)
## empty containers
fil2_list=list()
fil1_list=list()
ltlist=list()
gtlist=list()
eqlist=list()
####
  ############################################################ 
  ###data set t o use for x axis and y axis, correlation between x and y values will be used for dot size


  if(on_x == 'gee'){ x_data_set=geem}
  if(on_y == 'gee'){ y_data_set=geem}

  if(on_x == 'gei'){ x_data_set=geim}
  if(on_y == 'gei'){ y_data_set=geim}

  if(on_x == 'exd'){ x_data_set=dist_exom}
  if(on_y == 'exd'){ y_data_set=dist_exom}

  if(on_x == 'exd'){ x_data_set=dist_intm}
  if(on_y == 'exd'){ y_data_set=dist_intm}
  ###
  ### subset data to be used for plotting according to the list created after filtering above.
if(length(fil_list)>0){
    plot_x_data=x_data_set[fil_list,]
    plot_y_data=y_data_set[fil_list,]
}else{
   plot_x_data=c(NA,NA,NA)
   plot_y_data=c(NA,NA,NA)
}
 
## calculate correlation coefficient between all the exons and introns and put it in a data.frame s3d
library(Hmisc) ## for rcorr function to calculate corr coefficient with data containing NA
#tcor<-rcorr(as.numeric(plot_x_data[,2]),as.numeric(plot_y_data[,2]))


s3d<-data.frame()
for(ci in 1:Num_introns){
  for( cj in 1:Num_introns){
    
    #s3d<-rbind(s3d,c(ci,cj,cor(as.numeric(geemgt60[,ci]),as.numeric(geimgt60[,cj]))))
    ## rcorr{Hmisc} can calculate correlation from data with NA. cor does not do that
    s3d<-rbind(s3d,c(ci,cj,rcorr(plot_x_data[,ci],plot_y_data[,cj])$r[1,2]))
  }
}
## add column names to the data frame for the axis to display. 
colnames(s3d)<-c(labx,laby,"Corr coeff")
  
## create the title and subtitles to print below the plots.
  #subt<-paste("Coeff(",pars[div_col],"<",Num_introns,"):",  round(ltcor,2),"    Coeff(",pars[div_col],">",Num_introns,"):",  round(gtcor,2),sep=" ")
subt=paste("Num genes#",length(fil_list))
  maint<- substitute(expr = italic(m),env = list(m=name_list[i]))
  
### plot the data s3d by column1 and two onaxis and assign the column three (correlation) value to the size of dots.
  plot(s3d[,1],s3d[,2],xlab=colnames(s3d)[1],ylab=colnames(s3d)[2],cex=abs(s3d[,3]*10),pch=20,col = ifelse(s3d[,3] >= coeff_cutoff,'red',ifelse(s3d[,3] >= 0,'green',ifelse(s3d[,3] <= (-1*coeff_cutoff),'blue','yellow'))),main=maint,sub=subt)
  
}
dev.off()

print(paste("Saved image as",outfile,sep=" "))

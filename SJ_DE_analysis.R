## perform Differential gene analysis on junctions
overlaps<-function(rangex1,rangex2,rangey1,rangey2){
    frst = sort(c(as.integer(rangex1),as.integer(rangex2)))
    scnd = sort(c(as.integer(rangey1),as.integer(rangey2)))
    return(min(frst[2],scnd[2]) - max(frst[1],scnd[1]) > 0)
}

find_overlapp<-function(x,depth){
 #x="contig_1515:21836:21921:evm.TU.contig_1515.4:evm.TU.contig_1515.4 "
 #depth=5
  x=sort(x)
  u_ov_list=list()
  #print(paste("Recieved list of length",length(x)))
  for (i in 1:length(x)){
    for (d in 1:abs(depth)){
      if(i == d ) next
      if((i - d) < 0) next
      if(x[i] == x[i-d]) next
      #print(paste("Processing x: at i:",i,"at d:",d,x[i],sep=" "))
      junc1=unlist(strsplit(x[i]   ,':'))
      junc2=unlist(strsplit(x[i-d] ,':'))
      if(!(junc1[1] == junc2[1])) next
      #print(paste("Checking overlap between:",junc1[2],junc1[3],junc2[2],junc2[3],sep=" "))
      if(overlaps(junc1[2],junc1[3],junc2[2],junc2[3])){
        #print (paste("found overlap:",x[i],x[i-d], sep =" "))
        u_ov_list=c(u_ov_list,x[i],x[i-d])
      }
    }
  }
  
  u_ov_list=unique(unlist(u_ov_list,use.names = FALSE))
  return(u_ov_list)
}


group_deseq<-function(expdata,htsize,tiss_cdata_all,grp1,grp2){
  
  print(paste("Analyzing DE between groups:",grp1," and ",grp2))
  
  ### empty all variables to avoid accidental mixup
  tiss_data=tiss_cdata=tiss_dds=tiss_ddsFilt=tisssize=tisssamp=tiss_deseq_res=tiss_res_0.05 =tiss_res_padj0.05=tiss_ddsFilt_count=""
  
  ## Create new datasets
  tiss_cdata<-tiss_cdata_all[grep(grp1,tiss_cdata_all[,2]),]
  tiss_cdata<-rbind(tiss_cdata,tiss_cdata_all[grep(grp2,tiss_cdata_all[,2]),])
  
  ### create datasets based on tiss_cdata and perform analysis
  tiss_data<-as.matrix(expdata[tiss_cdata[,1]])
  
  colnames(tiss_cdata)<-c("Sample","tissue")
  tiss_dds<-DESeqDataSetFromMatrix(countData = tiss_data, colData = data.frame(tiss_cdata), design = ~ tissue)
  colnames(tiss_dds) <- colnames(tiss_data)
  tiss_ddsFilt<-tiss_dds[rowSums(counts(tiss_dds)) > 1,]
  tiss_ddsFilt<-estimateSizeFactors(tiss_ddsFilt)
  
  tisssamp<-colData(tiss_ddsFilt)$Sample
  
  ## use library size factors from htseq data.
  sizeFactors(tiss_ddsFilt)<-htsize[colnames(tiss_ddsFilt)]
  ##############################################################################################################
  ### all  
  tisssize<-sizeFactors(tiss_ddsFilt)
  tiss_ddsFilt<-DESeq(tiss_ddsFilt)
  tiss_deseq_res<-results(tiss_ddsFilt)
  tiss_ddsFilt_count=counts(tiss_ddsFilt,normalized=T)
  tiss_deseq_res_Alldata=cbind(as.matrix(tiss_deseq_res[rownames(tiss_deseq_res),]),as.matrix(tiss_ddsFilt_count[rownames(tiss_deseq_res),]))
  ofile=paste(grp1,"vs",grp2,"Diamond_DESEq_result.all.table",sep="_")
  write.table(tiss_deseq_res_Alldata,  file=ofile)
  print(paste("Saved results All DE in file:", ofile,sep=" "))
  ##############################################################################################################
  ### padj <=0.05
  tiss_res_padj0.05 <- tiss_deseq_res[! is.na(tiss_deseq_res$padj),]
  tiss_res_padj0.05 <- tiss_res_padj0.05[tiss_res_padj0.05$padj <= 0.05,]
  tiss_res_padj0.05_Alldata=cbind(as.matrix(tiss_res_padj0.05[rownames(tiss_res_padj0.05),]),as.matrix(tiss_ddsFilt_count[rownames(tiss_res_padj0.05),]))
  ofile=paste(grp1,"vs",grp2,"Diamond_DESEq_result.padj0.05.table",sep="_")
  write.table(tiss_res_padj0.05_Alldata,  file=ofile)
  #write.table(tiss_ddsFilt_count[rownames(tiss_res_padj0.05),],  file=paste(grp1,"vs",grp2,"Diamond_DESEq_normcount.padj0.05.table",sep="_"))
  print(paste("Saved results of padj > 0.05 in file:", ofile,sep=" "))
  ##############################################################################################################
  #### pvalue <=0.05  
  tiss_res_pval0.05 <- tiss_deseq_res[! is.na(tiss_deseq_res$pvalue),]
  tiss_res_pval0.05 <- tiss_res_pval0.05[tiss_res_pval0.05$pvalue <= 0.05,]
  tiss_res_pval0.05_Alldata=cbind(as.matrix(tiss_res_pval0.05[rownames(tiss_res_pval0.05),]),as.matrix(tiss_ddsFilt_count[rownames(tiss_res_pval0.05),]))
  ofile=paste(grp1,"vs",grp2,"Diamond_DESEq_result.pval0.05.table",sep="_")
  write.table(tiss_res_pval0.05_Alldata,file=ofile)
  print(paste("Saved results of pval > 0.05 in file:", ofile,sep=" "))
  
  
}




## test ing purpose only
#find_overlapp(swd_res_0.05_col,5)
##########################################



library(DESeq2)
wd="/home/ratnesh.singh/Turf/Turf_RNASeq_alternatesplicing"
res_fd=paste(wd,"DESeq_run_MergedStringtieGTF",sep="/")
setwd(wd)
sj_file = "Merged.SJ.againstStringtieGFF.UniqMappedRead.AllSplicetypes.table"

#### Using htseq count get real library size extimate
htseq_file = "/home/ratnesh.singh/Turf/htseq_results_all_mergedStringtieGTF/Turf_all_htseq_counts_all.csv"
htdata<-read.csv(htseq_file, header=T)
colnames(htdata)<-sub(".htseq","",(sub("STAR_map_","",colnames(htdata),perl=T,ignore.case = T)))
htd<-data.frame(htdata[,-1],row.names=htdata[,1])
htd<-htd[,order(colnames(htd))]
htcol<-colnames(htd)

##create aritficial groups to imitate real data to get library size.
htcold<-data.frame(cbind(htcol,c(rep("untreated",47),rep("treated",47))))
colnames(htcold)<-c("Sample","treatment")
htdds<-DESeqDataSetFromMatrix(countData = htd,colData = htcold, design = ~treatment)
colnames(htdds)<-colnames(htd)
## get the size factor for each library.
htdds<-estimateSizeFactors(htdds)
htsize<-sizeFactors(htdds)
htsamp<-colData(htdds)$Sample

##### Use SJ data and artificial groups to get normalized data for further analysis.
sjdata=read.table(sj_file,header=T)
colnames(sjdata)<-sub("_L000__R1R2_PSJ.out.tab","",sub("STAR_map_","",colnames(sjdata)))
expdata=data.frame(sjdata[,-1:-5],row.names=paste(sjdata[,1],sjdata[,2],sjdata[,3],sjdata[,4],sjdata[,5],sep=":"))
expdata<-expdata[,order(colnames(expdata))]
sjcol<-colnames(expdata)
cdata=cbind(sjcol,c(rep("salt",47),rep("water", 47)))
colnames(cdata)<-c("Sample","treatment")
ddsMat<-DESeqDataSetFromMatrix(countData = expdata, colData = data.frame(cdata), design = ~ treatment)
ddsMatFilt<-ddsMat[rowSums(counts(ddsMat)) > 1,]
ddsMatFilt<-estimateSizeFactors(ddsMatFilt)
sjsize<-sizeFactors(ddsMatFilt)
sjsamp<-colData(ddsMatFilt)$Sample

####### combine htsize and sjsize in one matrix
bsize<-cbind('htsize'=htsize[sjsamp],'sjsize'=sjsize[sjsamp])
rownames(bsize)<-sjsamp
#show(ddsMatFilt)
#sjdds<-DESeq(ddsMatFilt)
#htdds<-DES(htdds)

sj_fpm<-fpm(ddsMatFilt)
ht_fpm<-fpm(htdds)
## calculate SJ fpm based on library size determined by HT size.
sj_dds_HTsize<-ddsMatFilt
sizeFactors(sj_dds_HTsize)<-htsize

sj_fpm_HTsize<-fpm(sj_dds_HTsize)
colnames(sj_fpm_HTsize)<-colData(sj_dds_HTsize)$Sample
#####################################################################################################
### run DE analysis Pairwise ########################################################################
## Salt vs Water in Diamond
setwd(res_fd)
sjcol<-colnames(expdata)
sal_wat<-sjcol[grep("(Sal|Wat)",sjcol)]
swd_col<- sal_wat[grep("Diamond",sal_wat)]

swd_data<-expdata[swd_col]
swd_des<-c("Salt","Salt","Wat","Wat")

swd_cdata=cbind(swd_col,swd_des)
colnames(swd_cdata)<-c("Sample","treatment")
swd_dds<-DESeqDataSetFromMatrix(countData = swd_data, colData = data.frame(swd_cdata), design = ~ treatment)

swd_ddsFilt<-swd_dds[rowSums(counts(swd_dds)) > 1,]
##use htsize for analysis
swd_ddsFilt<-estimateSizeFactors(swd_ddsFilt)
sizeFactors(swd_ddsFilt)<-htsize[colnames(swd_ddsFilt)]
swdsize<-sizeFactors(swd_ddsFilt)



swdsamp           <-colData(swd_ddsFilt)$Sample
swd_ddsFilt       <-DESeq(swd_ddsFilt)
swd_deseq_res     <-results(swd_ddsFilt)
swd_res_0.05      <-swd_deseq_res[ ! is.na(swd_deseq_res$padj),]
swd_res_0.05      <-swd_deseq_res[ swd_res_0.05$padj <=0.05,]
swd_ddsFilt_count  =counts(swd_ddsFilt,normalized=T)
swd_res_0.05_row   =rownames(swd_res_0.05)
swd_res_0.05_OJlist=find_overlapp(swd_res_0.05_row,5)
swd_res_0.05_OJExp =as.matrix(swd_res_0.05[swd_res_0.05_OJlist,])
swd_res_0.05_OJCnt =as.matrix(swd_ddsFilt_count[swd_res_0.05_OJlist,])
swd_res_0.05_OJdata=cbind(swd_res_0.05_OJExp,swd_res_0.05_OJCnt)

swd_res_0.05_Alldata=cbind(as.matrix(swd_res_0.05[rownames(swd_res_0.05),]),as.matrix(swd_ddsFilt_count[rownames(swd_res_0.05),]))
            

write.table(swd_res_0.05_Alldata,file="Salt_water_Diamond_DESEq_result.table")

#write.table(swd_res_0.05_Alldata,                   file="Salt_water_Diamond_DESEq_result_normcount.padj0.05.table")
#write.table(swd_res_0.05_OJdata,  file="alt_water_Diamond_DESEq_result_normcount.OJunc.padj0.05.table")




#######################################################################################################################
### run DE analysis for pairwise in tissues
sjcol<-colnames(expdata)
tiss<-sjcol[grep("(Bud|Callus|Flower|Root|Rhizome|Diamond|Plant|Shoot|X)",sjcol)]
tiss_col<- tiss


### create tissue type file description
tiss_des=tiss_col
tiss_des=sub("Diamond_Rhizo\\.|Rhizome_\\.", "Rhizome_",tiss_des)
tiss_des=sub("Root\\.", "Root_",tiss_des)
tiss_des=gsub("Sal_Diamond|Wat_Diamond|X1|X2", "leaf_",tiss_des)
tiss_des=gsub("\\.", "_",tiss_des)
tiss_des=gsub("_\\w*", "",tiss_des)
tiss_cdata_all=cbind(tiss_col,tiss_des)
tiss_cdata=tiss_cdata_all
grp<-tiss_des[!duplicated(tiss_des)]

### create pairwise combination of tissues and run DE seq
i=1
j=2
for (i in 1:length(grp)){
  for (j in i:length(grp)){
    if(i == j)next

    group_deseq(expdata,htsize,tiss_cdata,grp[i],grp[j])
    
  }
}


### 
### create tissue type file description
trt_des=tiss_col
trt_des=sub("Diamond_Rhizo[\\d\\.]*Control_", "Rhizome.Control_",trt_des,perl=TRUE)
trt_des=sub("Diamond_Rhizo[\\d\\.C]*_[ATGC]*", "Rhizome.Cold",trt_des,perl =TRUE)
trt_des=sub("Root\\.", "Root_",trt_des)
trt_des=gsub("Sal_Diamond_[\\w\\d\\.]*", "Leaf.Salt_",trt_des)
trt_des=gsub("Wat_Diamond_[\\w\\d\\.]*", "Leaf.Water_",trt_des)
trt_des=gsub("_[\\w\\_\\.\\d]*", "",trt_des,perl=TRUE)
trt_des

trt_cdata_all=cbind(tiss_col,trt_des)
trt_cdata=trt_cdata_all
tgrp<-trt_des[!duplicated(trt_des)]
tgrp<-tgrp[grep("Leaf|Rhizome.C",tgrp)]

group_deseq(expdata,htsize,trt_cdata,"Rhizome.Control","Rhizome.Cold")
group_deseq(expdata,htsize,trt_cdata,"Leaf.Water","Leaf.Salt")





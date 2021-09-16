## Rscript: ctDNA_monitoring_stats.R
## Author: Michael Heskett
## Editing Author: Kami Chiotti [02.04.19]
## Usage: Rscript <path/to/infile.txt> <path/to/outdir/> <optional.cuffoff>
 
### NOTE: This script requires the 'dplyr' and 'sfsmisc' libraries to be installed.
# If using R locally, you may need to provide the location of your library, by
# uncommenting and editing the following example line:
#lib="/home/users/chiotti/R/x86_64-redhat-linux-gnu-library/3.5"
# input file must have columns in the following order: sample_id chromosome  position neg_ctrl_wt neg_ctrl_mut sample_wt sample_mut tissue_of_origin
 
#args=commandArgs(TRUE)
 
 
infile="infile.txt"
outdir="outdir"
 
if (is.na(cutoff)) {
  cutoff=0.005
}
 
library(dplyr)
library(sfsmisc)
 
### FUNCTIONS ###
 
get_overlap_coef=function(cmut, cwt, smut, swt, id){  ##version a/o 190129
  cmu=as.numeric(cmut/(cmut+cwt))
  smu=as.numeric(smut/(smut+swt))
  cvar=(as.numeric(cmut)*as.numeric(cwt))/(((cmut+cwt)^2)*(cmut+cwt+1))
  svar=(as.numeric(smut)*as.numeric(swt))/(((smut+swt)^2)*(smut+swt+1))
 
  if (max(c(cmut,smut)) == 0) {
    result=NA
    return(result)
  } else {
    xs=seq(max(0,min(cmu - 4*cvar, smu - 4*svar)), # xs must stay positive and go between 0,1
           min(1,max(cmu + 4*cvar, smu + 4*svar)),
           length.out=10000) ## this fails if mean and variance are 0 for both sample and case
  }
 
  f1=as.numeric(dbeta(xs, shape1=cmut, shape2=cwt))
  f2=as.numeric(dbeta(xs, shape1=smut, shape2=swt))
  minimum_curve=pmin(f1,f2)
  rslt=integrate.xy(xs,minimum_curve,min(xs),max(xs))
  if (max(f1,f2)==Inf) {
    ymax=2000
  } else {
    ymax=max(f1,f2)
  }
 
  sample_name=id
  jpeg(paste(outdir,id,'.jpeg',sep=""))
  plot(xs, f1, type="l", ylim=c(0, ymax), xlim=c(0, min(max(xs),1)) ,ylab="density")
  lines(xs, f2, lty="dotted")
  title(sample_name)
  dev.off()
  return(rslt)
}
### END OF FUNCTIONS ###
 
#### Per site
dat=read.table(infile,sep="\t",header=TRUE)
cn=colnames(dat)
colnames(dat)=c("id","chr","pos","cmut","cwt","smut","swt","tissue")
pv=('')
fil=('')
cvaf=('')
svaf=('')
for (i in 1:nrow(dat)){
  sample_id=paste(dat$id[i],"_",dat$pos[i],sep='')
  if(dat$cmut[i]==0) {
    dat$cmut[i]=1
  }
  cvaf[i]=(dat$cmut[i]/dat$cwt[i])*100
  svaf[i]=(dat$smut[i]/dat$swt[i])*100
  if(cvaf[i] > cutoff){
    fil[i]="FILTER"
  } else {
    fil[i]="KEEP"
  }
  pv[i]=get_overlap_coef(dat$cmut[i],dat$cwt[i],dat$smut[i],dat$swt[i],sample_id)
}
 
dat=cbind(dat,round(as.numeric(cvaf),4),round(as.numeric(svaf),4),round(as.numeric(pv),4),fil)
out=subset(dat,dat$fil=="KEEP")[,-ncol(dat)]
colnames(out)=c(cn,"negative_control_VAF","sample_VAF","raw_pvalue")
write.table(out,paste(outdir,unlist(strsplit(sample_id,"_"))[1],"_overlap_stats.txt",sep=''),sep="\t",row.names=F,col.names=T,quote=F)
 
agg.all=aggregate(out[,4:7], by=list(out$sample_id,rep("all",nrow(out))), sum)
agg.tis=aggregate(out[,4:7], by=list(out$sample_id,out$tissue_of_origin), sum)
agg=rbind(agg.all,agg.tis)
colnames(agg)=c("sample_id","tissue_type",colnames(agg[,3:ncol(agg)]))
 
raw_pvalue=('')
fil=('')
negative_control_VAF=('')
sample_VAF=('')
for (i in 1:nrow(agg)){
  sample_id=paste(agg$sample_id[i],"_",agg$tissue_type[i],sep='')
  negative_control_VAF[i]=(agg$neg_ctrl_mut[i]/agg$neg_ctrl_wt[i])*100
  sample_VAF[i]=(agg$sample_mut[i]/agg$sample_wt[i])*100
  raw_pvalue[i]=get_overlap_coef(agg$neg_ctrl_mut[i],agg$neg_ctrl_wt[i],agg$sample_mut[i],agg$sample_wt[i],sample_id)
}
 
 
tout=cbind(agg,round(as.numeric(negative_control_VAF),4),round(as.numeric(sample_VAF),4),round(as.numeric(raw_pvalue),4))
#tout=cbind(agg,round(as.numeric(negative_control_VAF),4),round(as.numeric(sample_VAF),4),round(as.numeric(raw_pvalue),4),fil)
#out=subset(dat,dat$fil=="KEEP")[,-ncol(dat)]
colnames(tout)=c(colnames(agg),"negative_control_VAF","sample_VAF","raw_pvalue")
write.table(tout,paste(outdir,unlist(strsplit(sample_id,"_"))[1],"_overlap_per_tissue_stats.txt",sep=''),sep="\t",row.names=F,col.names=T,quote=F)
 
q()
 
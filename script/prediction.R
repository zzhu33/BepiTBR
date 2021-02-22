##########  set up environment  ############

library("glmnet")

args = commandArgs(trailingOnly=TRUE)
#script_path="/project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/script"
#output_path="/project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/example/test_output_BepiTBR"
script_path=args[1]
output_path=args[2]
  
models=list(
  "bepipred1.0"=list(sd=0.427893141331944,
    model="cutoff1.5_bundary180_bin20_mid2mid_bepipred1.0_MixMHC2pred_seed449_alpha0.8_lambda0.021_modelsAndData.RData"),
  "bepipred2.0"=list(sd=0.394521694466676,
    model="cutoff1.5_bundary180_bin20_mid2mid_bepipred2.0_MixMHC2pred_seed449_alpha1_lambda0.021_modelsAndData.RData"),
  "LBEEP"=list(sd=0.599619428100351,
    model="cutoff2_bundary180_bin20_mid2mid_LBEEP_MixMHC2pred_seed449_alpha0.7_lambda0.001_modelsAndData.RData")
)

#######  prepare data  #######################

bin1.5=read.csv(list.files(paste(output_path,"/binned_csv1.5",sep=""),
  pattern="\\.csv",full.names = T),stringsAsFactors = F,row.names = 1)
bin2=read.csv(list.files(paste(output_path,"/binned_csv2",sep=""),
  pattern="\\.csv",full.names = T),stringsAsFactors = F,row.names = 1)

predictions=matrix(NA,ncol=2*length(models)+1,nrow=dim(bin1.5)[1])
colnames(predictions)=c(paste("base",names(models),sep="_"),
  paste("enhanced",names(models),sep="_"),"ensemble")
predictions=as.data.frame(predictions)
rownames(predictions)=rownames(bin2)

#########  predictions  #################

tepi="MixMHC2pred"

for (bepi in names(models))
{
  if (bepi=="LBEEP") {bin=bin2} else {bin=bin1.5}
  bin=bin[,!grepl("D_",colnames(bin))]
  bin=bin[,c(bepi,colnames(bin)[grepl(tepi,colnames(bin))])]
  bin=as.matrix(cbind(bin,bin[,1]*bin[,-1]))
  
  # base model
  predictions[,paste("base",bepi,sep="_")]=bin[,1]
  
  # enhanced model
  load(paste(script_path,"/",models[[bepi]]$model,sep=""))
  predictions[,paste("enhanced",bepi,sep="_")]=
    predict(fitEnhanced,newx=bin)[,1]
}

# ensemble
predictions$ensemble=1.8+rowMeans(sapply(names(models),function(x)
  predictions[,paste("enhanced",x,sep="_")]/models[[x]]$sd))

#########  write output  ##########

# polish output format
predictions$epitope=NA
for (file in list.files(paste(output_path,"/Bepi",sep="")))
{
  motif=read.table(paste(output_path,"/Bepi/",file,"/motif.fasta",sep=""),
    stringsAsFactors = F)
  predictions[file,"epitope"]=motif[2,1]
}

keep=colnames(predictions)!="epitope"
predictions[,keep]=round(predictions[,keep],d=4)

# write to file
file=paste(output_path,"/predictions.txt",sep="")
unlink(file)

if (all(predictions$base_bepipred2.0==0))
  {stop("BepiPred 2.0 must have failed!")}
if (all(predictions$base_bepipred1.0==0))
  {stop("BepiPred 1.0 must have failed!")}

write.table(predictions,file=file,
  row.names = T,col.names = T,quote=F,sep="\t")


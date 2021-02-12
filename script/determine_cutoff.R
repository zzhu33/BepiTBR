# run BepiTBR predictions using Ind-positive.txt 
# (the ind_validation_subset.txt subset) under
# /project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/script
# put results into /project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/example/test_output_BepiTBR

# read data
library("cutpointr")
pred=read.table("/project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/example/test_output_BepiTBR/predictions.txt",
  stringsAsFactors = F)

subset=read.table("/project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/script/ind_validation_subset.txt",
  stringsAsFactors = F)
pred=pred[rownames(pred) %in% subset[,1],]

# calculate ROC and cutoff
pred$label=1*grepl("Positive",rownames(pred))
cp=cutpointr(pred,ensemble,label,method=maximize_metric, metric=youden)

summary(cp)
pdf(file="/project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/example/test_output_BepiTBR/ROC.pdf",
  height=4,width=5.5)
plot(cp)
plot_metric(cp)
dev.off()

# best cutoff is -0.8 (rounded)
table(pred$ensemble>-0.8)
table(pred$label)

# base_bepipred1.0 as control
cp=cutpointr(pred,base_bepipred1.0,label,method=maximize_metric, metric=youden)
pdf(file="/project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/example/test_output_BepiTBR/ROC_bp1.pdf",
    height=4,width=5.5)
plot(cp)
dev.off()

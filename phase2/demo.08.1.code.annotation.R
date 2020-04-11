# -*- coding: utf-8 -*-
# @File    : demo.08.1.code.annotation.R
# @License : Copyright(C), 3DMed Co., Ltd.  
# @Author  : wenzhuan.xie@3dmedcare.com
# @Time    : 2020/04/04 13:42
# @IDE     : RStudio
# @Desc    : 需根据运行结果，对代码进行注释；
# @Desc    ：总结本代码目的：


#### 缓存清理 ####
rm(list = ls())
gc()
options(stringsAsFactors = F)

#### R包导入 ####
# 若报R包未安装，请根据phase1案例进行安装
library(magrittr)
library(pipeR)
library(ggpubr)
library(dplyr)
# library(clusterProfiler)
# library(org.Hs.eg.db)
library(biomaRt)
library(survminer)
library(survival)
library(openxlsx)

#### 变量定义 ####
cus.gene = 'FTO'
cutoff.method = 'median'

#### 载入数据 ####
tcga.clinical <- read.xlsx("./phase2/demo.08.file.cdr.xlsx",sheet = 1,rowNames = T)

tcga.clinical <- tcga.clinical %>%
  dplyr::rename(Tumor_Sample_Barcode = bcr_patient_barcode) # Anno:

table(tcga.clinical$type)

load("./phase2/demo.08.file.rnaseq.normal.Rdata") 
dim(rnaExpr)
rnaExpr[1:5,1:5]

clin.prad <- tcga.clinical %>%
  dplyr::filter(type == 'LUAD') # Anno:

#### 数据预处理 ####
dat.Expr <- as.data.frame(rnaExpr)
dat.Expr$Ensembl_ID <- rownames(dat.Expr) # Anno:

#' @1.ID转换
# 定义好数据集
listMarts()
cus.mart <- useMart('ENSEMBL_MART_ENSEMBL')
all.sets <- listDatasets(cus.mart) 
cus.sets <- useDataset('hsapiens_gene_ensembl',mart = cus.mart)
# 进行转换
id.att <- listAttributes(cus.sets)
id.deg <- getBM(attributes = c('ensembl_gene_id','chromosome_name','hgnc_symbol','entrezgene_id','gene_biotype'), # 需需要转换未何种格式
                filters = 'ensembl_gene_id',    # 指定数据的输入类型
                values = dat.Expr$Ensembl_ID,   # Anno:
                mart = cus.sets)

id.df <- data.frame(table(id.deg$gene_biotype))

id.selet <- id.deg %>%  
  dplyr::filter(gene_biotype %in% c('protein_coding','lncRNA'), # Anno:
                hgnc_symbol != '') %>% 
  dplyr::distinct(hgnc_symbol,.keep_all = T) %>% 
  dplyr::rename(Ensembl_ID = ensembl_gene_id)

#' @2.数据合并
dat.plot <- merge(id.selet,dat.Expr,by='Ensembl_ID')# Anno:
dat.plot <- dat.plot[,-c(1,2,4,5)]                  # Anno:
dat.plot.mat <- dat.plot[,-c(1)]
dat.plot$median <- rowMedians(as.matrix(dat.plot.mat))   # Anno:
dat.plot <- arrange(dat.plot,desc(median))               # Anno:
dat.plot <- distinct(dat.plot,hgnc_symbol,.keep_all = T) # Anno:
dat.plot <- dplyr::select(dat.plot,-median)              # Anno:
rownames(dat.plot) <- dat.plot[,1]
dat.plot <- dat.plot[,-1]

sur.df <- as.data.frame(t(dat.plot[c(cus.gene),]))       # Anno:
sur.df$barcode <- substr(rownames(sur.df),1,12)          # Anno:
sur.df$type <- substr(rownames(sur.df),14,15)            # Anno:
sur.df <- sur.df[which(sur.df$type == '01'),]            # Anno:

cln.df <- dplyr::select(clin.prad,barcode = Tumor_Sample_Barcode,times = OS.time,status = OS) # Anno:
cln.df <- cln.df[which(cln.df$times >= 30),] # Anno:
cln.df$times <- cln.df$times / 30            # Anno:

sur.df <- merge(cln.df,sur.df,by='barcode')



#' @计算cutoff
if ( cutoff.method == 'median') { 
  best.cutoff = median(sur.df$FTO)     # Anno:
}else{
  res.cut <- surv_cutpoint(sur.df, 
                           time = "times",  
                           event = "status",
                           variables = "FTO")
  best.cutoff <- res.cut$cutpoint[[1]] # Anno:
  cutoff.method = 'optimal'
}

#' @分组
sur.df$risk <- ifelse(sur.df$FTO >= best.cutoff,'High','Low') # Anno:
sur.df$risk <- factor(sur.df$risk,levels = c('Low','High'))   # Anno:
table(sur.df$status)
#' @保存
write.table(sur.df,file = paste0("./phase2/demo.08.Out.",cus.gene[1],".",cutoff.method,".1.tsv"),
            sep = "\t",col.names = T,row.names = F)



#' @kmplot
plot.km = function(plotdata){
  #' @plotdata：为数据框，必含几下几个字段：times（时间维度为Month）,status, risk(分类值，如High和low)
  {
    data.survdiff <- survdiff(Surv(plotdata$times, plotdata$status, type = 'right') ~ plotdata$risk, data = plotdata)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    x = summary(coxph(Surv(times, status,type = 'right')~risk, data = plotdata,method = 'efron'))
    HR = signif(x$coef[2], digits=2)
    up95 = signif(x$conf.int[,"upper .95"],2)
    low95 = signif(x$conf.int[,"lower .95"], 2)
    HR <- paste("HR:", round(HR,2), sep = "")
    CI <- paste("95%CI: ", paste(round(low95,2), round(up95,2), sep = "-"), sep = "")
    HRCI <- paste0(HR,"(",CI,")")
  }
  sfit <- survfit(Surv(times,status) ~ risk,data = plotdata)  
  p <- ggsurvplot(sfit, 
                  data = plotdata,
                  conf.int=F, #置信区间
                  pval.coord = c(max(plotdata$times) * 0.35, 0.8), # 位值信息，可调整
                  pval=paste(HRCI,pval = ifelse(p.val < 0.001, "p < 0.001", paste("p = ",round(p.val,4), sep = "")), sep = "\n"),
                  palette = "jama",
                  risk.table =T, 
                  ncensor.plot = F,
                  legend.labs = c(paste0(toupper(names(table(plotdata$risk)[1])),"(",table(plotdata$risk)[[1]],")"),
                                  paste0(toupper(names(table(plotdata$risk)[2])),"(",table(plotdata$risk)[[2]],")")))+  # 根据factor指定的顺序给定改字符串顺序
    guides(color=guide_legend(override.aes=list(fill=NA))) + 
    labs(x = paste0("Months"),y='Percent survival')
  p$table <- p$table + theme(axis.line = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title.y = element_blank(),
                             axis.title.x = element_blank(),
                             axis.text.x = element_blank(),
                             plot.title = element_text(size=12)) 
  return(p)
}
p <- plot.km(sur.df)

pdf(file = paste0("./phase2/demo.08.Out.",cus.gene[1],".",cutoff.method,".1.pdf"),
    onefile = F)
print(p)
dev.off()

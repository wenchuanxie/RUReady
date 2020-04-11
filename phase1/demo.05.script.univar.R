#' @File    : demo.05.script.univar.R
#' @License : Copyright(C), 3DMed
#' @Author  : wenzhuan.xie@3dmedcare.com
#' @Time    : 2020/02/15
#' @IDE     : RStudio
#' @Desc    : 寻找连续变量的最优分组cutoff值

rm(list = ls())
gc()
options(stringsAsFactors = F)

if (!require(survival)) BiocManager::install('survival',update = F)
if (!require(dplyr)) BiocManager::install('dplyr',update = F)

library(survival)
library(dplyr)

#### 数据准备 ####
plot.df <- read.table("./phase1/demo.05.file.univar.tsv",
                      sep = "\t",
                      header = T)  # 包含行名、times和status以及待分析得基因；基因表达值经过TMM标准化（也可直接log2）

#### 单因素分析 ####
train <- plot.df
Sur <- Surv(train$times, train$status)                 # 建立生存对象
uni.cox <- function(x){
  fml <- as.formula(paste0('Sur~', x))
  gcox <- coxph(fml, train)
  cox_sum <- summary(gcox)
  HR <- round(cox_sum$coefficients[,2],2)
  PValue <- round(cox_sum$coefficients[,5],6)
  CI <- paste0(round(cox_sum$conf.int[,3:4],2),collapse='-')
  uni.res <- data.frame('Characteristics' = x,
                        'Hazard Ratio' = HR,
                        'CI95' = CI,
                        'P value' = PValue)
  return(uni.res)
}

uni.list <- lapply(setdiff(colnames(train), c("times","status")),uni.cox)
uni.df <- plyr::ldply(uni.list, data.frame)
uni.sig <- na.omit(uni.df[uni.df$P.value <= 0.05,]) # 筛选单因素分析结果小于等于设定阈值的变量

write.table(uni.sig,file = "./phase1/demo.05.Out.univar.tsv",sep = "\t",col.names = T,row.names = F,quote = F)




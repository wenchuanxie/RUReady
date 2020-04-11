#' @File    : demo.02.script.cutoff.R
#' @License : Copyright(C), 3DMed
#' @Author  : wenzhuan.xie@3dmedcare.com
#' @Time    : 2020/02/15
#' @IDE     : RStudio
#' @Desc    : 寻找连续变量的最优分组cutoff值

rm(list = ls())
gc()
options(stringsAsFactors = F)

if (!require(survminer)) BiocManager::install('survminer',update = F)
if (!require(dplyr)) BiocManager::install('dplyr',update = F)

library(survminer)
library(dplyr)

#### 数据准备 ####
plot.df <- read.table("./phase1/demo.02.file.survival.tsv",
                      sep = "\t",
                      header = T)

#### 寻找时间依赖的最优cutoff ####
res.cut <- surv_cutpoint(plot.df, 
                         time = "OS_MONTHS",  # 随访时间
                         event = "OS_STATUS", # 值为0（事件未发生）或者1（事件发生）
                         variables = c("Immune",'Stromal','ESTIMATE')) # 需要分组的连续变量，多个以,分隔
summary(res.cut) # 查看各个变量的最优分界值

# 绘图查看最佳分组cutoff
plot(res.cut, "Immune", palette = "npg") # Plot cutpoint
plot(res.cut, "Stromal", palette = "npg") # Plot cutpoint
plot(res.cut, "ESTIMATE", palette = "npg") # Plot cutpoint
# 根据最佳cutoff定义分组
res.cat <- surv_categorize(res.cut)  
head(res.cat)

# 将定义分组的数据与原始数据合并
plot.df.2 <- plot.df %>%
  dplyr::select(PATIENT_ID) %>%
  dplyr::bind_cols(res.cat)

write.table(plot.df.2,file = "./phase1/demo.02.Out.cutoff.tsv",sep = "\t",col.names = T,row.names = F,quote = F)



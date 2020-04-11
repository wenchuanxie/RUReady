#' @File    : demo.01.script.Table1.R
#' @License : Copyright(C), 3DMed
#' @Author  : wenzhuan.xie@3dmedcare.com
#' @Time    : 2020/02/15
#' @IDE     : RStudio
#' @Desc    : 创建临床三线表

rm(list = ls())
gc()
options(stringsAsFactors = F)

if (!require(BiocManager)) install.packages("BiocManager")
library(BiocManager)

if (!require(arsenal)) BiocManager::install('arsenal',update = F)
if (!require(dplyr)) BiocManager::install('dplyr',update = F)

library(arsenal)
library(dplyr)

#### 数据准备 ####
load("./phase1/demo.01.file.clinical.Rdata")

#### 绘制Table1 ####
table_one <- tableby(Type ~ AGE + SEX + RACE + WEIGHT + TUMOR_STAGE + RADIATION_THERAPY, 
                     data = plot.df)                # ~ 前的属性表示分组变量，~ 后的属性为需要汇总的变量
summary(table_one)                                  # 在Console窗口打印显示
write2word(table_one, "./phase1/demo.01.Out.clinical.doc") # 输出为word文件
write2(table_one,"./phase1/demo.01.Out.clinical.xls")     # 输出为excel文件


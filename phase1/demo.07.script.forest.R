#' @File    : demo.07.script.forest.R
#' @License : Copyright(C), 3DMed
#' @Author  : wenzhuan.xie@3dmedcare.com
#' @Time    : 2020/02/15
#' @IDE     : RStudio
#' @Desc    : 多因素分析森林图

rm(list = ls())
gc()
options(stringsAsFactors = F)


if (!require(survival)) BiocManager::install('survival',update = F)
if (!require(survminer)) BiocManager::install('survminer',update = F)
if (!require(ggplot2)) BiocManager::install('ggplot2',update = F)
if (!require(patchwork)) BiocManager::install('patchwork',update = F)

library(survival)
library(survminer)
library(ggplot2)
library(patchwork)
#### 载入数据 ####
dat.df <- read.table("./phase1/demo.07.file.forest.tsv",
                     sep = "\t",header = T)
 
#### 多因素cox回归 ####
# Surv：用于创建生存数据对象
# coxph：构建COX回归模型
# 1.1 对三个因素：sex、rx、adhere进行回归分析
class(dat.df$sex) 
table(dat.df$sex)
class(dat.df$rx)
table(dat.df$rx)
class(dat.df$adhere)
table(dat.df$adhere)  # 三个都是分类变量
model <- coxph( Surv(time, status) ~ sex + rx + adhere, data = colon )
summary(model)

# 1.2 绘制森林图
p1 <- ggforest(model) 
print(p1)


# 2.1 对三个因素：sex、rx、adhere进行回归分析
dat.df <- within(dat.df, {
  sex <- factor(sex, labels = c("female", "male"))
  rx <- factor(rx,levels = c("Obs","Lev","Lev+5FU"))
  adhere <- as.character(adhere)
  differ <- factor(differ, labels = c("well", "moderate", "poor"))
  extent <- factor(extent, labels = c("submuc.", "muscle", "serosa", "contig."))
})
model2 <- coxph(Surv(time, status) ~ sex + rx + adhere + differ + extent,data = dat.df )
summary(model2)

# 2.2 绘制森林图
p2 <- ggforest(model2) 
print(p2)

# 3 合并两图，对比有发现差别吗？
p <- (p1 / p2) + plot_layout(ncol = 1,heights = c(1,2)) + plot_annotation(tag_levels = c("A","B"))
print(p)

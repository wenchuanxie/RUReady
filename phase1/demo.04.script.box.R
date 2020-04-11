#' @File    : demo.04.script.box.R
#' @License : Copyright(C), 3DMed
#' @Author  : wenzhuan.xie@3dmedcare.com
#' @Time    : 2020/02/15
#' @IDE     : RStudio
#' @Desc    : 寻找连续变量的最优分组cutoff值

rm(list = ls())
gc()
options(stringsAsFactors = F)

if (!require(ggpubr)) BiocManager::install('ggpubr',update = F)
if (!require(dplyr)) BiocManager::install('dplyr',update = F)

library(ggpubr)
library(dplyr)

#### 数据准备 ####
plot.df <- read.table("./phase1/demo.04.file.box.tsv",
                      sep = "\t",
                      header = T)
#### 根据分组绘制两组箱线图 ####
plot.score <- plot.df %>% 
  tidyr::drop_na(TUMOR_STAGE)  # 剔除TUMOR_STAGE值为NA得记录
p <- ggboxplot(plot.score, x = "TUMOR_STAGE",
               y = c('Stromal','Immune','ESTIMATE'),
               combine = TRUE,
               color = "TUMOR_STAGE", palette = "jama",
               ylab = "Score", 
               xlab = "",
               add = "jitter",                              # Add jittered points
               add.params = list(size = 1, jitter = 0.3)    # Point size and the amount of jittering
              ) + 
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) +
  stat_compare_means(method = "anova")  # 三组以上时，使用方差分析
print(p)
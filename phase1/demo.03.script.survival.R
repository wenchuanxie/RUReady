#' @File    : demo.03.script.survival.R
#' @License : Copyright(C), 3DMed
#' @Author  : wenzhuan.xie@3dmedcare.com
#' @Time    : 2020/02/15
#' @IDE     : RStudio
#' @Desc    : 寻找连续变量的最优分组cutoff值

rm(list = ls())
gc()
options(stringsAsFactors = F)

if (!require(survminer)) BiocManager::install('survminer',update = F)
if (!require(survival)) BiocManager::install('survival',update = F)
if (!require(dplyr)) BiocManager::install('dplyr',update = F)

library(survminer)
library(survival)
library(dplyr)

#### 数据准备 ####
plot.df <- read.table("./phase1/demo.02.Out.cutoff.tsv",
                      sep = "\t",
                      header = T)
#### 根据分组绘制两组生存曲线 ####
# 构建绘图函数
plot.km = function(plotdata){
  #' @plotdata：为数据框，必含几下几个字段：times（时间维度为Month）,status, risk(分类值，如High和low)
  {
    # R coxph中，两组（risk：mut vs wt）比较，默认首字母在前的组是ref组（mut），另一组是对比组（wt），即后者相对于前者比较，是什么结果（HR的解释）
    # 若需要反过来比较，即mut比wt是什么结果，则执行以下语句（或者通过factor指定，本处在调用函数前设定好ref，因此不执行下述语句）：
    # plotdata$risk  <- plotdata$risk %>% 
    #   as.factor() %>% 
    #   forcats::fct_rev()
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
                  pval.coord = c(max(plotdata$times) * 0.02, 0.1), # 位值信息，可调整
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

# 绘制Immune评分分组的生存曲线
km.df <- plot.df %>%
  dplyr::rename(times = OS_MONTHS,status = OS_STATUS,
                risk = Immune) %>%
  dplyr::mutate(risk = factor(risk,levels = c('high','low'))) # ***因子化分组,前者(high)为ref
p <- plot.km(km.df)
print(p)
gg <- plotly::ggplotly(p$plot)
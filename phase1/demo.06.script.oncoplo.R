# @File    : customOncoplot.R
# @License : Copyright(C), 3DMed
# @Author  : wenzhuan.xie@3dmedcare.com
# @Time    : 2019/08/20 08:42
# @IDE     : RStudio
# @Desc    : 瀑布图绘制


# 安装包
if(!require("BiocManager")){
  install.packages("BiocManager")
}
if(!require("dplyr")){
  BiocManager::install("dplyr",update = F,ask = F)
}  
if(!require("data.table")){
  BiocManager::install("data.table",update = F,ask = F)
} 
if(!require("magrittr")){
  BiocManager::install("magrittr",update = F,ask = F)
}
if(!require("GenVisR")){
  BiocManager::install("GenVisR",update = F,ask = F)
}
if(!require("maftools")){
  BiocManager::install("maftools",update = F,ask = F)
}

# 导入包
library(dplyr)
library(GenVisR)
library(magrittr)
library(data.table)
library(maftools)

# 清空已有的对象
rm(list = ls())  # Attention：本操作会清空所有的数据
gc()
options(stringsAsFactors = F)


########################################################################
## 瀑布图制作方法一
## 要求事先准备好复核要求的数据文件
# 文件内容格式：
# (1)包含三列，列名分别为：sample, gene, variant_class
# (2)sample：样本的ID，或者其他能标识样本的符号（非中文）
# (3)gene: 变异的基因
# (4)variant_class: 基因变异类型
# 以当前demo.06.file.oncoplot-1.tsv文件为例(tsv后缀的文件可以用excel打开)，前3行数据：
#      sample       gene        variant_class
#  1807187807       BRAF                indel
#  1807187807     MAP3K1       non-synonymous
#  1807187807      SLIT2       non-synonymous

# !!! 注意文件不要保存为xls或者xlsx格式 !!!

# 载入硬盘中存放的文件数据
data <- data.table::fread("./phase1/demo.06.file.oncoplot-1.tsv",encoding = "UTF-8")  # 读入一个数据文件
dim(data)                                                   # 查看数据维度，如果列数比行数多很多，则不执行下一句View指令
#View(data)                                                  # 查看数据内容
# 修正每一列的属性
data <- data %>%
  dplyr::mutate(sample = as.character(sample),
                gene = as.character(gene),
                variant_class = as.character(variant_class))

# 绘制突变瀑布图
vars <- as.character(unique(data$variant_class))
waterfall(data,
          #maxGenes = 10,      # 表示绘制丰度最高的前10个基因
          plotGenes = c("RET","SPEN","APC","STK11","ALK","ROS1","BRCA2"),# 表示绘制指定的基因，与上述 maxGenes 二选一
          fileType = "Custom",
          variant_class_order = vars)



########################################################################
## 瀑布图制作方法二
## 要求事先准备好复核要求的数据文件
# 文件内容格式：
# (1)包含5列，列名分别为：sample, gene, c_dot,p_dot,variant_type
# (2)sample：样本的ID，或者其他能标识样本的符号（非中文）
# (3)gene: 变异的基因
# (4)c_dot: 碱基变化
# (5)p_dot: 氨基酸变化
# (6)variant_type: 基因变异分级
# 以当前demo.06.file.oncoplot-2.tsv文件为例(tsv后缀的文件可以用excel打开)，前3行数据：
#   sample      gene       c_dot       p_dot     variant_type
#  6180S01      KRAS      c.G34T      p.G12C              SNV
#  6180S01     NTRK2     c.C568A     p.Q190K              SNV
#  6180S01       HGF    c.C1361A     p.T454K              SNV

# !!! 注意文件不要保存为xls或者xlsx格式 !!!

# 载入硬盘中存放的文件数据
data <- data.table::fread("./phase1/demo.06.file.oncoplot-2.tsv",encoding = "UTF-8")
dim(data)                                                   # 查看数据维度，如果列数比行数多很多，则不执行下一句View指令
#View(data)                                                  # 查看数据内容

# 根据 variant_type、c_dot、p_dot判断variant_class
input.data <- data
input.data$variant_class <- ifelse(input.data$variant_type == "SNV" | input.data$variant_type == "Germline",
                               ifelse(grepl("fs",input.data$c_dot) | grepl("fs",input.data$p_dot),
                                      "frameshift indel",
                                      ifelse(grepl("ins",input.data$c_dot) | grepl("del",input.data$c_dot) | 
                                               grepl("ins",input.data$p_dot) | grepl("del",input.data$p_dot) | 
                                               grepl("dup",input.data$p_dot),
                                             "non-frameshift indel",
                                             ifelse(grepl("ext",input.data$p_dot),
                                                    "nonstop",
                                                    ifelse(grepl("^-",substring(input.data$c_dot,3)),
                                                           "5UTR",
                                                           ifelse(grepl("\\+",input.data$c_dot) | grepl("\\-",input.data$c_dot),
                                                                  "splice",
                                                                  ifelse(grepl("\\*$",input.data$p_dot),
                                                                         "nonsense",
                                                                         ifelse(substr(input.data$p_dot,3,3) == substr(input.data$p_dot,nchar(input.data$p_dot),nchar(input.data$p_dot)),
                                                                                "synonymous",    # p.G123G, 表示同义突变
                                                                                "missense") # 否则为非同义突变
                                                                  )
                                                           )
                                                    )
                                             )
                                      )
                               ),
                               ifelse(input.data$variant_type == "CNV",
                                      ifelse(input.data$p_dot == "gain",
                                             "gain",
                                             "loss"
                                      ),
                                      ifelse(input.data$variant_type == "fusion",
                                             "fusion",
                                             "others")
                               )
                               )
data <- input.data 

# 修正每一列的属性
data <- data %>%
  dplyr::mutate(sample = as.character(sample),
                gene = as.character(gene),
                variant_class = as.character(variant_class))

# 绘制突变瀑布图
vars <- as.character(unique(data$variant_class))
waterfall(data,
          maxGenes = 10,      # 表示绘制丰度最高的前10个基因
          #plotGenes = c("EGFR","RET","SPEN","APC","STK11","ALK","ROS1","BRCA2"),# 表示绘制指定的基因，与上述 maxGenes 二选一
          fileType = "Custom",
          variant_class_order = vars)


########################################################################
## 瀑布图制作方法三：使用原始变异名称，自定义颜色显示
## 要求事先准备好复合要求的数据文件:demo.06.file.oncoplot-3.tsv
## 强制字段：Hugo_Symbol、Chromosome、Start_Position、End_Position、Reference_Allele、Tumor_Seq_Allele2、Variant_Classification、Variant_Type以及Tumor_Sample_Barcode。
## 可选但建议包含的字段：VAF（Variant Allele Frequency）以及氨基酸变化信息。
if(T){
  rm(list = ls())
  gc()
  # 读取基因变异数据
  raw.data <- read.table("./phase1/demo.06.file.oncoplot-3.tsv",sep = "\t",header = T,stringsAsFactors = F)
  # 根据 variant_type、c_dot、p_dot判断variant_class
  input.data <- raw.data
  input.data$Variant_Classification <- ifelse(input.data$Variant_Type == "SNV" | input.data$Variant_Type == "Germline",
                                     ifelse(grepl("fs",input.data$c_dot) | grepl("fs",input.data$p_dot),
                                            "frameshift indel",
                                            ifelse(grepl("ins",input.data$c_dot) | grepl("del",input.data$c_dot) | 
                                                     grepl("ins",input.data$p_dot) | grepl("del",input.data$p_dot) | 
                                                     grepl("dup",input.data$p_dot),
                                                   "non-frameshift indel",
                                                   ifelse(grepl("ext",input.data$p_dot),
                                                          "nonstop",
                                                          ifelse(grepl("^-",substring(input.data$c_dot,3)),
                                                                 "5UTR",
                                                                 ifelse(grepl("\\+",input.data$c_dot) | grepl("\\-",input.data$c_dot),
                                                                        "splice",
                                                                        ifelse(grepl("\\*$",input.data$p_dot),
                                                                               "nonsense",
                                                                               ifelse(substr(input.data$p_dot,3,3) == substr(input.data$p_dot,nchar(input.data$p_dot),nchar(input.data$p_dot)),
                                                                                      "synonymous",    # p.G123G, 表示同义突变
                                                                                      "missense") # 否则为非同义突变
                                                                        )
                                                                 )
                                                          )
                                                   )
                                            )
                                     ),
                                     ifelse(input.data$Variant_Type == "CNV",
                                            ifelse(input.data$p_dot == "gain",
                                                   "gain",
                                                   "loss"
                                            ),
                                            ifelse(input.data$Variant_Type == "fusion",
                                                   "fusion",
                                                   "others")
                                     )
  )
  raw.data <- input.data 
  raw.data <- raw.data %>%
    dplyr::distinct(Tumor_Sample_Barcode,Hugo_Symbol,.keep_all = T) # 按样本名和基因名随机去重复。慎重操作
  
  write.table(raw.data,file = "./demo.06.file.oncoplot-3.maf",quote = F,sep = "\t",col.names = T,row.names = F)
  
  library(maftools)
  vc_nonSyn <- setdiff(as.character(unique(raw.data$Variant_Classification)),'synonymous')
  raw.maf = read.maf(maf = "./phase1/demo.06.file.oncoplot-3.maf",
                     vc_nonSyn = vc_nonSyn  # 自定义非同义突变名称，因为上述类型不是MAF指定的类型名称
                     )
  #oncoplot for top ten mutated genes.
  {
    colnames(raw.maf@variant.classification.summary)
    # 自定义颜色
    library(RColorBrewer)
    display.brewer.all()
    #One can use any colors, here in this example color palette from RColorBrewer package is used
    vc_cols = RColorBrewer::brewer.pal(n = length(vc_nonSyn), name = 'Paired')
    names(vc_cols) = c(vc_nonSyn)
    print(vc_cols)
  }
  library(devEMF)
  emf(file = "./phase1/demo.06.Out.oncoplot.emf",emfPlus = F,width =9 ,height = 6,family = 'Arial')
  oncoplot(maf = raw.maf, 
           top = 10, # 与参数'genes'二选一
           #genes = c('EGFR','KRAS','ROS1','STK11','ATM','RET','MET'),
           #showTumorSampleBarcodes = T, # 显示样本名
           colors = vc_cols)
  dev.off()
  
}

















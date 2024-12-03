library(microeco);
library(ggsci);
library(aplot);
library(tidyverse);
library(reshape2);
library(ggpubr);
library(viridis)

source("代码作者微信cgxr410/cgxr410.R")
# 示例数据格式

OTU_sample_data1=otu_table_16S  ;  
Species_sample_data1=taxonomy_table_16S

OTU_sample_data2=otu_table_ITS  ; 
Species_sample_data2=taxonomy_table_ITS


#—————————————↓类型一：数据符合示例数据的格式,直接读取↓—————————————————————————————————————————————————————
otu <- read.csv('otu_table.csv', row.names = 1);
taxonomy <- read.csv('taxonomy.csv', row.names = 1)

#—————————————↓类型二：不符合示例数据格式,进行数据处理，获得示例数据的格式↓————————————————————————————————————————————
data <- read.csv('otu_table_taxonomy.csv', row.names = 1)
data1 <- VX_cgxr410_Data_preprocessing(VXcgxr410)
otu=data1$otu ; taxonomy=data1$taxonomy
#————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
#设定分组
group_labels <- factor(c(rep("A", 4),rep("AA", 4), rep("B", 4))) 
#FAPROTAX功能预测并作图
num1=30 ; num2=20 ; width1 = 6 ; height1 = 4 ; size.word=6  #图1参数
y.size=7;x.size=7; width2 = 15 ; height2 = 8 ;#图2参数
method = "anova" #可选项：2组分析用 t.test、wilcox.test；3组用 anova、kruskal.test、chisq.test


color = c("#4DBBD5FF", "#E64B35FF", "#00A087FF", "#3C5488FF") #可根据分组数调整,颜色数可大于分组数，不能小于分组数

#执行分析1:FAPROTAX原核微生物功预测
result_FAPROTAX <- VX_cgxr410_FAPROTAX(VXcgxr410)
#执行分析2:FUNGuild真核微生物功预测
result_FUNGuild <- VX_cgxr410_FUNGuild(VXcgxr410)
result_FungalTraits <- VX_cgxr410_FungalTraits(VXcgxr410)






setwd("E:\\桌面\\microeco包妙用之网络鲁棒性、易损性、内聚力和网络稳定性")


🌟基础操作之meconetcomp包等依赖包安装、数据微表创建、共现性网络构建
#############包安装和加载###################
# 安装所需的 R 包
# aplot: microeco 包中 trans_venn 类的依赖包
# agricolae: 用于 Duncan 的新多重范围检验
packages <- c("meconetcomp", "rgexf", "pheatmap", "aplot", "agricolae")
# 检查每个包是否已安装，如果未安装则安装
for(x in packages) {
  if(!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE) # 自动安装包及其所有依赖项
  }
}
#安装并加载ape包
if(!require("ape")){install.packages("ape")}
library(ape)
library(microeco)
library(meconetcomp)
library(magrittr)
library(igraph)
library(ggplot2)
library(rgexf)
library(ggpubr)
library(agricolae)
theme_set(theme_bw())
data(soil_amp)
##################数据基本操作#######################
# 加载示例数据；16S rRNA基因扩增子测序数据
#读取样本分组信息表
sample_info_16S <- read.csv("sample_table.csv", row.names = 1)
#读取特征表
otu_table_16S <- read.csv("feature_table.csv",row.names = 1)
#读取系统发育树
phylo_tree_16S <- read.tree("phylo_tree.tre", tree.names = NULL)
#读取物种分类表
taxonomy_table_16S <- read.csv("tax_table.csv",row.names = 1)
# 让我们创建需要信息的 microtable 对象
soil_amp <- microtable$new(sample_table = sample_info_16S, otu_table = otu_table_16S, tax_table = taxonomy_table_16S, phylo_tree = phylo_tree_16S)

##微生物共现性网络构建
#########创建三组网络############
# 选择"IW"组的样本
# 使用clone对soil_amp对象深度拷贝
tmp <- clone(soil_amp)
# 直接更改样本表
tmp$sample_table %<>% subset(Group == "IW")
# 清理对象中的所有文件
tmp$tidy_dataset()
# 使用filter_thres参数过滤低丰度特征
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
# 设置 p 值阈值
# 设置相关系数阈值
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
# 将网络加入列表中
soil_amp_network$IW <- tmp


# 选择"TW"组的样本
tmp <- clone(soil_amp)
tmp$sample_table %<>% subset(Group == "TW")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
soil_amp_network$TW <- tmp


# 选择"CW"组的样本
tmp <- clone(soil_amp)
tmp$sample_table %<>% subset(Group == "CW")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0005)
tmp$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
soil_amp_network$CW <- tmp
# 现在已创建网络列表 soil_amp_network




##网络鲁棒性指标计算
################## 网络鲁棒性 ##################
# 创建一个 'robustness' 对象来评估网络的鲁棒性，即网络在不同移除策略下的抗扰动能力
# remove_strategy 参数指定移除策略，包括随机移除边("edge_rand")、移除最强的边("edge_strong")、随机移除节点("node_rand")和移除度数高的节点("node_degree_high")
# remove_ratio 指定移除比例，从0到0.99，步长为0.1
# measure 参数指定鲁棒性指标，Eff 表示效率，Eigen 表示特征值，Pcr 表示聚类系数
# run = 10 表示每种配置运行10次，便于计算平均鲁棒性
tmp <- robustness$new(soil_amp_network, remove_strategy = c("edge_rand", "edge_strong", "node_rand", "node_degree_high"), 
                      remove_ratio = seq(0, 0.99, 0.1), measure = c("Eff", "Eigen", "Pcr"), run = 10)
# 查看鲁棒性分析的结果表，包括每次实验的详细输出
View(tmp$res_table)
write.csv(tmp$res_table,"robustness_detail.csv")
# 查看鲁棒性分析的总结，包括每种策略的鲁棒性统计结果
View(tmp$res_summary)
write.csv(tmp$res_summary,"robustness_summary.csv")
# 绘制鲁棒性分析图表，显示不同移除策略下网络鲁棒性随移除比例的变化
g1<- tmp$plot(linewidth = 1)
g1
ggsave("robustness.pdf", g1, width = 10, height = 7)





##网络易损性指标计算
################节点易损性 #########################
# 使用 'vulnerability' 函数计算每个节点的易损性，即节点被移除时对网络整体稳定性的影响
# 易损性指标有助于识别关键节点
vul_table <- vulnerability(soil_amp_network)
View(vul_table)  # 查看易损性结果表，包含每个节点的易损性评分
#输出所有的节点的易损性值
write.csv(vul_table,"vulnerability_all_otu.csv")




##网络内聚力及稳定性
############### 内聚力 ####################
# 创建 'cohesionclass' 对象来计算网络的凝聚力，即网络中节点之间的紧密程度
t1 <- cohesionclass$new(soil_amp_network)
# 查看凝聚力计算的结果，包括样本级别的凝聚力信息
View(t1$res_list$sample)
#输出内聚力计算结果
write.csv(t1$res_list$sample,"Cohesion.csv")
# 查看连通性信息，例如基于各特征的连接性分析结果
View(t1$res_list$feature)
#输出连通性计算结果
write.csv(t1$res_list$feature,"连通性.csv")
# 计算不同处理组间凝聚力的差异，使用 ANOVA 检验
t1$cal_diff(method = "anova")
# 绘制内聚力图表，使用正连通性（r_pos）作为指标，显示各组别的凝聚力水平
g1<-t1$plot(measure = "r_pos")
g1
ggsave("r_pos.pdf", g1, width = 7, height = 7)


#基于正负内聚力的网络稳定性计算
net<- t1$res_list$sample
library(dplyr)
# 创建新列，值为 c_neg 的绝对值除以 c_pos
net <- net %>%
  mutate(stability = abs(c_neg) / c_pos)  # 计算并添加新列
write.csv(net,"network stability.csv")


























































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
if(!require("ape")){install.packages("ape")}

setwd("E:\\桌面\\microeco包-网络属性（拓扑、节点）及多网络比较（子网络拓扑）")

###########数据读取和构建微表#############
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


# 首先创建一个列表
soil_amp_network <- list()

#本示例提供了3个组别，可根据自己的数据进行设计。
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




#网络模块化计算
soil_amp_network %<>% cal_module(undirected_method = "cluster_fast_greedy")

#网络拓扑属性计算
tmp <- cal_network_attr(soil_amp_network)
# tmp 是一个包含网络属性的数据框,输出tmp文件
write.csv(tmp,"tmp2.csv")

#提取所有网络的节点和边属性
soil_amp_network %<>% get_node_table(node_roles = TRUE) %>% get_edge_table


####多网络特性比较####
#跨网络比较节点
#比较网络重叠节点和差异节点
# 使用 meconetcomp 包中的 node_comp 函数将所有网络中的节点转换为一个新表对象，然后使用 trans_venn 类轻松分析节点重叠情况。
tmp <- node_comp(soil_amp_network, property = "name") # 获取节点分布
tmp1 <- trans_venn$new(tmp, ratio = "numratio") # 获取节点交集
g1 <- tmp1$plot_venn(fill_color = FALSE) # 绘制韦恩图
g1
ggsave("soil_amp_node_overlap.pdf", g1, width = 7, height = 6)

# 计算 Jaccard 距离以反映网络的整体差异
tmp$cal_betadiv(method = "jaccard")
tmp$beta_diversity$jaccard


####跨网络比较边####
# 获取跨网络的边分布
tmp <- edge_comp(soil_amp_network)
tmp1 <- trans_venn$new(tmp, ratio = "numratio") # 获取边的交集
g2 <- tmp1$plot_venn(fill_color = FALSE)
g2
ggsave("soil_amp_edge_overlap.pdf", g2, width = 7, height = 6)

# 计算 Jaccard 距离
tmp$cal_betadiv(method = "jaccard")
tmp$beta_diversity$jaccard


####将网络重叠的边提取到新网络####
# 获取边的分布和交集
tmp <- edge_comp(soil_amp_network)
tmp1 <- trans_venn$new(tmp)
# 将交集结果转换为 microtable 对象
tmp2 <- tmp1$trans_comm()
# 提取所有三组网络 ("IW", "TW" 和 "CW") 的交集
Intersec_all <- subset_network(soil_amp_network, venn = tmp2, name = "IW&TW")
Intersec_all$save_network("Intersec_all.gexf") # 保存为 gexf 格式

####比较边中成对节点的系统发育距离####
# 过滤低丰度特征以加速计算
node_names <- unique(unlist(lapply(soil_amp_network, function(x) { colnames(x$data_abund) })))
filter_soil_amp <- microeco::clone(soil_amp)
filter_soil_amp$otu_table <- filter_soil_amp$otu_table[node_names, ]
filter_soil_amp$tidy_dataset()

# 获取系统发育距离矩阵
phylogenetic_distance_soil <- as.matrix(cophenetic(filter_soil_amp$phylo_tree))

# 使用正负标签进行计算
tmp <- edge_node_distance$new(network_list = soil_amp_network, dis_matrix = phylogenetic_distance_soil, label = c("+", "-"))
tmp$cal_diff(method = "anova")
# 可视化
g3 <- tmp$plot(add = "none", add_sig = TRUE, add_sig_text_size = 5) + ylab("Phylogenetic distance")
g3
ggsave("soil_amp_phylo_distance.pdf", g3, width = 7, height = 6)

# 显示具有至少 10 个节点和正边的不同模块
tmp <- edge_node_distance$new(network_list = soil_amp_network, dis_matrix = phylogenetic_distance_soil, 
                              label = "+", with_module = TRUE, module_thres = 10)
tmp$cal_diff(method = "anova")
g4 <- tmp$plot(add = "none", add_sig = TRUE, add_sig_text_size = 5) + ylab("Phylogenetic distance")
g4
ggsave("soil_amp_phylo_distance_modules.pdf", g4, width = 8, height = 6)


####比较跨网络边的节点来源####
soil_amp_network_edgetax <- edge_tax_comp(soil_amp_network, taxrank = "Phylum", label = "+", rel = TRUE)
# 过滤小数量特征
soil_amp_network_edgetax <- soil_amp_network_edgetax[apply(soil_amp_network_edgetax, 1, mean) > 0.01, ]

# 可视化
g5 <- pheatmap::pheatmap(soil_amp_network_edgetax, display_numbers = TRUE)
g5
ggsave("soil_amp_edge_tax_comp.pdf", g5, width = 7, height = 7)


####比较子网络的拓扑性质####
#该部分首先根据soil_amp数据集中每个样本所含的OTU提取soil_amp_network中每个网络的子网络，然后计算子网络的全局拓扑性质，所有操作都封装在meconetcomp包的subnet_property函数中。
# 计算所有子网络的全局属性
tmp <- subnet_property(soil_amp_network)
# 为相关性分析准备数据
# 将第二列（样本名）作为行名
rownames(tmp) <- tmp[, 2]
# 删除前两列（网络名称和样本名称）
tmp <- tmp[, -c(1:2)]
# 加载现成的环境因子和多样性数据表
data(soil_measure_diversity)
# 创建 trans_env 对象，添加环境因子和多样性数据
tmp1 <- trans_env$new(dataset = soil_amp, add_data = soil_measure_diversity)
# 计算相关性矩阵，按组分析，使用添加的丰度表 tmp，Spearman 方法
tmp1$cal_cor(use_data = "other", by_group = "Group", add_abund_table = tmp, cor_method = "spearman")
# 生成相关性热图
g6 <- tmp1$plot_cor()
g6
# 保存热图至 PDF 文件
ggsave("soil_amp_subnet_property.pdf", g6, width = 11, height = 5)




参考：https://mp.weixin.qq.com/s/lz5m2tajO1eYQsWlHCy2ew 并做了修正
#plspm 包来进行 偏最小二乘路径建模（PLS-PM），一种用于多变量因果关系建模和分析的技术。PLS-PM 可用于探索潜在变量之间的因果关系，尤其适用于较少样本和高维度的数据。下面是对代码的详细解释：
1. 加载和准备数据
#plspm 包是用于PLS路径建模的R包，首先需要加载。
library(plspm) #加载plspm包
#read.csv() 用于加载一个CSV文件，假设该文件包含了OTU（操作分类单元）数据。row.names = 1 表示把数据集的第一列作为行名。
test_otu<-read.csv("test_otu.csv", row.names = 1) #读取数据集，是一个CSV文件

2. 构建数据框（df）
df <- test_otu[1:50, c(5:12, 1:4)]
#df 是一个包含 50 个样本和 12 个观测变量的数据框。test_otu[1:50, c(5:12, 1:4)] 选择了前 50 行数据，并选择了第 5 到第 12 列以及第 1 到第 4 列的数据。

colnames(df) <- c("SOC", "TN", "TP", "TK", "Bacteria1", "Bacteria2", "Fungi1", "Fungi2", "N_gene1", "N_gene2", "N_proc1", "N_proc2")
#colnames(df) 为数据框命名，定义了 12 个观测变量的名称，如“SOC”表示土壤有机碳，"Bacteria1" 和 "Bacteria2" 表示细菌种群，等等。

3. 可选的标准化
# df = scale(df[, 1:12], center = TRUE, scale = TRUE)
#如果需要，可以对数据进行标准化（即均值为 0，标准差为 1），特别是当数据的量级差异较大时，标准化是常用的步骤。

4. 定义内模型（inner model）
#内模型（Inner model）：内模型定义了潜在变量（如 Soil, Bac_div, Fun_div, N_gene, N_proc）之间的因果关系。这些潜在变量在路径模型中通过矩阵表示。
Soil <- c(0, 0, 0, 0, 0)
Bac_div <- c(1, 0, 0, 0, 0)
Fun_div <- c(1, 0, 0, 0, 0)
N_gene <- c(0, 1, 1, 0, 0)
N_proc <- c(0, 1, 1, 1, 0)
df_path <- rbind(Soil, Bac_div, Fun_div, N_gene, N_proc)

#每一行对应一个潜在变量的依赖关系。例如：
#Soil <- c(0, 0, 0, 0, 0) 表示 Soil 不受其他潜在变量的影响（即它是一个外部变量）。
#Bac_div <- c(1, 0, 0, 0, 0) 表示 Bac_div 受 Soil 的影响。
#df_path 是一个下三角矩阵，表示潜在变量之间的关系。

5. 定义外模型（outer model）
#外模型（Outer model）：外模型定义了每个潜在变量与其对应的观测变量（也就是模型的测量部分）之间的关系。
df_blocks = list(1:2, 3:4, 5:8, 9:10, 11:12) # 映射观测变量到潜在变量
df_modes = rep("A", 5)
#df_blocks 列出了每个潜在变量对应的观测变量的列索引。例如：第一个潜在变量 Soil 对应的观测变量是 SOC 和 TN（列 1 和列 2）。
#df_modes：指定了每个潜在变量的模式类型，"A" 表示这些潜在变量是反射性的（reflective），即潜在变量通过观测变量来定义。

6. 执行PLS-PM分析
df_pls = plspm(df, df_path, df_blocks, modes = df_modes, boot.val = TRUE)
#plspm() 函数用于进行PLS路径建模，返回的对象 df_pls 包含了模型分析的结果。
#df：数据框，包含了观测变量的数据。
#df_path：内模型矩阵，表示潜在变量之间的关系。
#df_blocks：外模型矩阵，表示每个潜在变量与其对应的观测变量之间的关系。
#modes：潜在变量的模式，"A" 表示反射模型。
#boot.val = TRUE：表示进行自助法（Bootstrap）抽样，用于计算统计显著性。

7. 可视化与结果分析
plot(df_pls, what="loadings", arr.width = 0.1, show.values = TRUE, lcol = 'gray')
#plot(df_pls, what="loadings", ...)：绘制潜在变量与观测变量之间的相关性（loading）。通过这些加载值，查看潜在变量和观测变量之间的关系，通常加载值要大于0.7才有较强的关联性。

innerplot(df_pls, colpos = '#CA5023', colneg = '#457CC3', show.values = TRUE, lcol = 'gray20', box.lwd = 0)
#innerplot(df_pls, ...)：绘制路径模型的可视化图，显示潜在变量之间的路径系数和相关性。

df_pls$inner_model
#df_pls$inner_model 输出的是 路径模型 中各个潜在变量之间的路径系数估计值及其统计显著性。每个潜在变量的路径系数估计值（Estimate）表示该变量对另一个潜在变量的影响强度，而其他列则提供该估计值的标准误差、t值和p值，用于衡量该路径系数的统计显著性。

df_pls$gof
#df_pls$gof：查看拟合优度（Goodness of Fit，GOF）指标，GOF值大于0.7表示模型拟合较好。

summary(df_pls)
#summary(df_pls)：查看模型的详细结果，包括路径系数、显著性和拟合度等信息。
#通过AI或者PTT对路径系数和显著性进一步可视化

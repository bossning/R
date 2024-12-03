VX_cgxr410_FAPROTAX <- function(VXcgxr410) {
  if (!dir.exists("16s")) {
    dir.create("16s")
  }
  m_table <- microtable$new(otu_table = otu ,tax_table = taxonomy)
  t2 <- trans_func$new(m_table)
  t2$cal_spe_func(prok_database = "FAPROTAX")
  t2$cal_spe_func_perc(abundance_weighted = FALSE)
  result_FAPROTAX =t2$res_spe_func_perc
  if (ncol(result_FAPROTAX) < num1) {
    num1 <- ncol(result_FAPROTAX)
  }
  if (ncol(result_FAPROTAX) < num2) {
    num2 <- ncol(result_FAPROTAX)
  }
  # 绘图准备：计算总和，排序，并选择前30列
  top30_columns <- names(sort(colSums(result_FAPROTAX, na.rm = TRUE), decreasing = TRUE)[1:num1])
  selected_data <- result_FAPROTAX[, top30_columns]
  # 标准化数据，并将行名转为列 'sample'
  data_frame <- as.data.frame(scale(selected_data), row.names = rownames(result_FAPROTAX))
  data_frame$sample <- rownames(data_frame)
  # 从宽格式转为长格式，并进行绘图
  long_data <- melt(data_frame, id.vars = 'sample')
  p <- ggplot(long_data, aes(x = sample, y = variable)) +
    geom_point(aes(size = abs(value), color = value)) +
    scale_size_area(max_size = 2) +
    scale_color_gsea() +
    labs(x = "", y = "") +
    theme_classic(base_size = size.word) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  # 保存图像和最终结果
  ggsave('16s/1.plot_FAPROTAX.pdf', plot = p, width = width1, height = height1)
  write.csv(t2$res_spe_func_perc, file = "16s/result_FAPROTAX.csv", row.names = TRUE)
  
  # 绘图准备：计算总和，排序，并选择前30列
  top15_columns <- names(sort(colSums(result_FAPROTAX, na.rm = TRUE), decreasing = TRUE)[1:num2])
  selected_data2 <- result_FAPROTAX[, top15_columns]
  # 将result_FAPROTAX的第一列移除，并转换为数据框
  selected_data2 <- selected_data2[, -1]
  selected_data2 <- as.data.frame(selected_data2)
  # 对数值进行标准化
  result_FAPROTAX_scaled <- as.data.frame(scale(selected_data2))
  # 添加分组信息
  result_FAPROTAX_scaled$Group <- group_labels
  # 转换数据格式并添加分组信息
  df <- melt(result_FAPROTAX_scaled, id.vars = "Group", variable.name = "Function", value.name = "score")
  # 计算最大值并加上5%
  max_value <- max(df$score)
  label_y_position <- max_value + max_value * 0.1
  # 绘制箱线图
  p1 <- ggplot(data = df, aes(x = Function, y = score, fill = Group)) +
    geom_boxplot(width = 0.3, position = position_dodge(0.5), outlier.colour = NA) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_fill_manual(values = color) +  # 根据你的组数调整颜色
    stat_compare_means(aes(group = Group), method = "anova", label = "p.signif", hide.ns = TRUE, label.y = c(label_y_position)) +
    theme(axis.text.x = element_text(size = x.size, colour = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(size = y.size, angle = 0),
          axis.title.x = element_text(size = 5),
          axis.title.y = element_text(size = 7),
          legend.text = element_text(size = 7)) +  # 控制图例中文本字体大小
    xlab("")  # 去掉 x 轴标签
  
  # 打印图形
  ggsave('16s/2.boxplot_FAPROTAX.pdf', plot = p1, width = width2, height = height2)
  
  return(result_FAPROTAX=result_FAPROTAX)
}

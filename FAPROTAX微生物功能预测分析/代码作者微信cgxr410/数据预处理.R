VX_cgxr410_Data_preprocessing <- function(VXcgxr410) {

# 提取taxonomy列，并保持行名
taxonomy <- data$taxonomy
# 保留原来的行名到新的数据框
taxonomy_df <- as.data.frame(taxonomy, row.names = row.names(data))
# 根据第一列进行分列，分列符号为;
taxonomy <- data.frame(do.call(rbind, strsplit(as.character(taxonomy_df$taxonomy), ";")))
clean_taxonomy <- function(tax) {
  return(gsub("_unclassified.*", "", tax))
}
taxonomy <- as.data.frame(lapply(taxonomy, clean_taxonomy), stringsAsFactors = FALSE)
row.names(taxonomy) <- row.names(taxonomy_df)
colnames(taxonomy)[1:7] <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# 删除taxonomy列
otu <- data[, -which(names(data) == "taxonomy")]
return(list(otu=otu,taxonomy=taxonomy))
}

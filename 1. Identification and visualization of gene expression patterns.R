library(DESeq2)
library(dplyr)
library(tidyverse)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(gridExtra)
raw_counts <- read.csv("GSE152418_GeneLevel_Raw_data.csv")
colnames(raw_counts)[1]<-"genes"
raw_counts<-aggregate(.~genes,raw_counts,mean)%>%
  column_to_rownames(.,"genes")

raw_dim<-dimnames(raw_counts)
raw_counts<-apply(raw_counts,2,as.integer)
dimnames(raw_counts)<-raw_dim
raw_counts <-as.data.frame(raw_counts)

# And here's a corresponding sample_info dataframe:
sample_info <- read.csv("GSE152418_filtered_metadata.csv")%>%
  select(.,c(1,6))%>%
  column_to_rownames(.,"ID")
sample_info_BOX <- cbind(sample=rownames(sample_info),sample_info)
raw_counts <- raw_counts[,colnames(raw_counts)%in%rownames(sample_info)]

# Loading DESeq2 library
library(DESeq2)

# Creating DESeq object
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = sample_info, design = ~ severity)
dds$severity <- relevel(dds$severity, ref = "Healthy")
dds <- DESeq(dds)
vst_data <- vst(dds)
normalized_counts <- counts(dds, normalized=TRUE)

# Here, assuming that 'severity' column in your sample_info dataframe is a factor variable with four levels: 'healthy', 'light', 'neutral', and 'heavy'



res_M_vs_H <- results(dds, contrast = c("severity", "Moderate", "Healthy"))
res_S_vs_M <- results(dds, contrast = c("severity", "Severe", "Moderate"))
res_I_vs_S <- results(dds, contrast = c("severity", "ICU", "Severe"))

df <- data.frame(Gene = rownames(res_M_vs_H),
                 log2FC_M_vs_H = res_M_vs_H$log2FoldChange,
                 p_M_vs_H = res_M_vs_H$pvalue,
                 log2FC_S_vs_M = res_S_vs_M$log2FoldChange,
                 p_S_vs_M = res_S_vs_M$pvalue,
                 log2FC_I_vs_S = res_I_vs_S$log2FoldChange,
                 p_I_vs_S = res_I_vs_S$pvalue)

df <- df[!is.na(df$p_M_vs_H) & 
           !is.na(df$p_S_vs_M) & 
           !is.na(df$p_I_vs_S), ]
df<-df[df$p_M_vs_H<0.5&df$p_S_vs_M<0.5&df$p_I_vs_S<0.5,]

# Define the pattern conditions
# 1. Genes that are continuously upregulated
upregulated_genes <- df[df$log2FC_M_vs_H > 0 & df$log2FC_S_vs_M > 0 & df$log2FC_I_vs_S > 0, ]

# 2. Genes that are initially upregulated and then downregulated
up_down_genes <- df[df$log2FC_M_vs_H > 0 & df$log2FC_S_vs_M < 0 & df$log2FC_I_vs_S < 0, ]

# 3. Genes that are initially downregulated and then upregulated
down_up_genes <- df[df$log2FC_M_vs_H < 0 & df$log2FC_S_vs_M > 0 & df$log2FC_I_vs_S > 0, ]

# 4. Genes that are continuously downregulated
downregulated_genes <- df[df$log2FC_M_vs_H < 0 & df$log2FC_S_vs_M < 0 & df$log2FC_I_vs_S < 0, ]

# Let's assume you want to display the first 10 genes in a boxplot
selected_genes <- c("COL1A2","CD177","CD63",# up
                    "CD28","CLOCK","CXCR4",# down
                    "CCL2","CXCL11","CD38",# updown
                    "STIM1","CD36","CXCL16"# downup
                    ) 

# load("all_results.Rdata")

# Let's assume you want to display the first 10 genes in a boxplot
# selected_genes <- common_genes[1:24]

# Get the variance stabilized data for the selected genes
selected_vst_counts <- assay(vst_data)[selected_genes, ]

# Convert the matrix to a data frame and reshape it
long_format_df <- as.data.frame(t(selected_vst_counts)) %>%
  rownames_to_column(var = "sample") %>%
  gather(gene, expression, -sample)

# Merge with sample_info to get the severity for each sample
long_format_df <- merge(long_format_df, sample_info_BOX, by = "sample")
long_format_df$severity <- factor(long_format_df$severity,levels = c("Healthy","Moderate","Severe","ICU"))
long_format_df <- long_format_df %>%
  mutate(gene_model = case_when(
    gene %in% upregulated_genes$Gene ~ "Consistent_Up",
    gene %in% downregulated_genes$Gene ~ "Consistent_Down",
    gene %in% up_down_genes$Gene ~ "Up_then_Down",
    gene %in% down_up_genes$Gene ~ "Down_then_Up",
    TRUE ~ NA_character_  # 默认条件，当其他条件都不满足时
  ))
long_format_df$gene_model <- factor(long_format_df$gene_model,levels = c("Consistent_Up","Consistent_Down","Up_then_Down","Down_then_Up"))

# Plot with ggplot2
# title = "Gene expression across different severities",
selected_plot <- ggplot(long_format_df, aes(x = severity, y = expression, fill = severity)) +
  facet_wrap(~ gene_model + gene, scales = "free_y", nrow = 4) +
  geom_boxplot() +
  # facet_wrap(~ gene, scales = "free_y") +
  theme_minimal() +
  labs( x = "Severity", y = "Scaled expression level", fill = "Severity")+
  scale_fill_manual(values = c("#42B540FF","#00468BFF","#e28822","#ED0000FF"))+
  theme(text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1))
selected_plot
ggsave(filename = "model_gene_plot.pdf",plot = selected_plot,width =10,height = 8 )

# 分别为每一个gene_model画一个图
plots_list <- lapply(unique(long_format_df$gene_model), function(gene_model) {
  df_subset <- subset(long_format_df, gene_model == gene_model)
  
  ggplot(df_subset, aes(x = severity, y = expression, fill = severity)) +
    geom_boxplot() +
    facet_wrap(~ gene, scales = "free_y") +
    theme_minimal() +
    labs(title = paste("Gene expression across different severities: ", gene_model), x = "Severity", y = "Scaled expression level", fill = "Severity")+
    scale_fill_manual(values = c("#42B540FF","#00468BFF","#e28822","#ED0000FF"))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
})

# 使用ggarrange和marrangeGrob进行分页
arranged_plots <- ggarrange(plotlist = plots_list, ncol = 2, nrow = 2)
# marrangeGrob(arranged_plots, nrow = 2, ncol = 2, pages = "as_needed")

# save(list = ls(),file = "all_results.Rdata")
# 
# load("all_results.Rdata")

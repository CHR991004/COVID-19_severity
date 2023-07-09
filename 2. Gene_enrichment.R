load("all_results.Rdata")
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
up_genes <- data.frame(genes=upregulated_genes$Gene,type = rep("Consis_up",nrow(upregulated_genes)))
down_genes <- data.frame(genes=downregulated_genes$Gene,type = rep("Consis_down",nrow(downregulated_genes)))
updown_genes <- data.frame(genes=up_down_genes$Gene,type = rep("Up_then_down",nrow(up_down_genes)))
downup_genes <- data.frame(genes=down_up_genes$Gene,type = rep("Down_then_up",nrow(down_up_genes)))
all_genes <- rbind(up_genes,down_genes,updown_genes,downup_genes)
all_genes$entrez <- mget(all_genes$genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)%>%
  as.character()
all_genes<- all_genes[all_genes$entrez!="NA",]
enrich_KEGG_obj <- compareCluster(entrez~type,
                                  data=all_genes,
                                  fun="enrichKEGG",qvalueCutoff = 1,pvalueCutoff = 0.05)
enrich_KEGG_obj <- setReadable(enrich_KEGG_obj, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
enrich_KEGG_obj@compareClusterResult$type <- factor(enrich_KEGG_obj@compareClusterResult$type,
                                                     levels = c("Consis_up","Up_then_down","Consis_down","Down_then_up"))

enrich_KEGG_obj@compareClusterResult <- enrich_KEGG_obj@compareClusterResult[order(enrich_KEGG_obj@compareClusterResult$type),]
KEGG_selected <- read.table("KEGG_selected.txt",sep = "\t")[,1]
enrich_KEGG_obj@compareClusterResult <- enrich_KEGG_obj@compareClusterResult[enrich_KEGG_obj@compareClusterResult$Description%in%KEGG_selected,]

KEGG_dot_plot <- enrichplot::dotplot(enrich_KEGG_obj,
        x="type",showCategory = 40)+
  # facet_grid(~num)+ # width可设置条形图宽度
  xlab("") +
  theme_bw()+
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE))+
  scale_color_gradient2(low = "#470c5f", high = "#fde725", mid = "#1f958b", 
                       midpoint = 0.025,  space = "Lab",na.value = "white",guide = "colourbar") +
  theme(axis.title = element_text(size = 15,face = "bold",colour = "black"), # 可以设置标题字体大小、是否加粗、颜色
        panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"), # 外围加框，fill为填充色
        axis.ticks = element_line(linewidth = 1),
        plot.title = element_text(size = 30, face = "bold"),
        axis.text.x = element_text(angle=45,hjust=1,size = 12),
        axis.text.y = element_text(size=12))
write.table(enrich_KEGG_obj@compareClusterResult,file = "KEGG_results_all.xls",sep = "\t",quote = F,row.names = F,col.names = T)
ggsave(filename = "KEGG_dotplot.PDF",plot = KEGG_dot_plot,height = 10,width = 8)

enrich_GO_obj <- compareCluster(entrez~type,
                                data=na.omit(all_genes),
                                fun="enrichGO", OrgDb='org.Hs.eg.db',ont="ALL")

enrich_GO_obj@compareClusterResult$type <- factor(enrich_GO_obj@compareClusterResult$type,
                                                    levels = c("Consis_down","Consis_up","Down_then_up","Up_then_down"))
write.table(enrich_GO_obj@compareClusterResult,file = "GO_results_manual.xls",sep = "\t",quote = F,row.names = F,col.names = T)
enrich_GO_obj@compareClusterResult <- enrich_GO_obj@compareClusterResult[order(enrich_GO_obj@compareClusterResult$ONTOLOGY),]
GO_selected<- unique(read.table("GO_selected_manual.txt",sep = "\t",header = F,check.names = F,quote = "",comment.char = "")[,1])
# GO_selected<- c(GO_selected,enrich_GO_obj@compareClusterResult[enrich_GO_obj@compareClusterResult$type=="Consis_down","Description"])
enrich_GO_obj@compareClusterResult <- enrich_GO_obj@compareClusterResult[enrich_GO_obj@compareClusterResult$Description%in%GO_selected,]
write.table(x = enrich_GO_obj@compareClusterResult,"GO_selected_final2.xls",sep = "\t",row.names = F,col.names = T,quote = F)
GO_dot_plot <- dotplot(enrich_GO_obj,
        x="type",showCategory = 10)+
  facet_grid(~ONTOLOGY)+ # width可设置条形图宽度
  xlab("") +
  theme_bw()+
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE))+
  scale_color_gradient2(low = "#470c5f", high = "#fde725", mid = "#1f958b", 
                        midpoint = 0.02,  space = "Lab",na.value = "white",guide = "colourbar") +
  theme(axis.title = element_text(size = 15,face = "bold",colour = "black"), # 可以设置标题字体大小、是否加粗、颜色
        panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"), # 外围加框，fill为填充色
        axis.ticks = element_line(linewidth = 1),
        plot.title = element_text(size = 30, face = "bold"),
        axis.text.x = element_text(angle=45,hjust=1,size = 12),
        axis.text.y = element_text(size=12))
ggsave(filename = "GO_dotplot.PDF",plot = GO_dot_plot,height = 8,width = 10)

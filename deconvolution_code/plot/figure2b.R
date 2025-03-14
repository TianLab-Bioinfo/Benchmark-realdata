setwd('example_dir')  
source('./code/plot/function_figure2b.R')

scRNA_list <- readRDS('./scRNA_list.rds')
ER_files = sort(list.files('./COAD/TCGA/normal/'))
TNBC_files = sort(list.files('./COAD/TCGA/tumor/'))

ER_dir = './COAD/TCGA/normal/'
TNBC_dir = './COAD/TCGA/tumor/'

key_trend_dir <- "./scRNA_standard/median_total/COAD.rds"
key_trend <- readRDS(key_trend_dir)
row.names(key_trend) <- key_trend[,1]
df2=trend_plot(ER_files,ER_dir,TNBC_dir,key_trend,'TCGA-COAD')

plot_df <- as.matrix(df2)

plot_df_melt <- melt(plot_df)
colnames(plot_df_melt) <- c('celltype','method','trend') 
plot_df_melt[is.na(plot_df_melt)] <- 0


key_trend$Trend <- factor(key_trend$Trend , levels = c('+',  0,'-'))
celltype_order <- key_trend[order(key_trend$Trend,key_trend$PValue),'Celltype']

plot_df_melt$celltype <- factor(plot_df_melt$celltype ,levels = celltype_order)
method_order <- c('Ground_True','DWLS','music','CIBERSORT',  'bisque','ReCIDE','bayes')
plot_df_melt$method <- factor(plot_df_melt$method, levels = method_order)


plot_df_melt$trend2 <- plot_df_melt$trend
plot_df_melt$trend <- as.factor(plot_df_melt$trend)
p <- ggplot(plot_df_melt,aes(x=celltype,y=method,fill=trend))+ 
  geom_tile(aes(fill = trend), color = "#F5F5F5", size = 0.1)+
  coord_equal()+  
  scale_fill_manual(values = c("-1"="#C1DCFF", "1"="#FFD8E8", "0.5" = "white","-0.5" = "white","0" = "white","0.3" = "#91C79A")) +
  geom_text(aes(label = ifelse(abs(trend2 ) == 0.5, ifelse(trend2  > 0, "+", "-"), "")), color = 'black', size = 4) +  # 添加符号
  # scale_fill_gradient2(low="#C1DCFF", high="#FFD8E8", mid="white")+   #"#87CEEB""#B1D3FF"
  #TCGAm_text(aes(label=round(value,2)),color='black',size = 2)+
  theme(legend.position = 'none',axis.title.x = element_blank(),axis.title.y =element_blank(),axis.text.x =element_text(angle =315,hjust =0.5,vjust = 0.5))+
  guides(
    fill = guide_colorbar(
      barwidth = 1,        
      barheight = 3,        
      title.position = "top", 
      title.hjust = 0.5   
      
    ))+coord_fixed(ratio = 1)+
  annotate("rect", xmin = 0.5, xmax = length(unique(plot_df_melt$celltype)) + 0.5,
           ymin = 0.5, ymax = length(unique(plot_df_melt$method)) + 0.5,
           color = "black", fill = NA, size = 0.7)+
  geom_text(aes(label = ifelse(abs(trend2) == 0.5, ifelse(trend2 > 0, "+", "-"), "")), color = 'black', size = 4) 
p



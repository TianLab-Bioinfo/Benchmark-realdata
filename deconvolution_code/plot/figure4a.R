setwd('example_dir')
file_dir5 = './COAD_df.rds'
file_plot <- readRDS(file_dir5)
file_plot=file_plot[!file_plot$celltype %in% c('P-value','Correlation','RMSE'),]
file_plot$p_GEO <- as.numeric(file_plot$p_GEO)
file_plot$p_TCGA <- as.numeric(file_plot$p_TCGA)
file_plot[file_plot$direction_TCGA=='-','TCGA'] <- 0-file_plot[file_plot$direction_TCGA=='-','p_TCGA']
file_plot[file_plot$direction_TCGA=='+','TCGA'] <- file_plot[file_plot$direction_TCGA=='+','p_TCGA']
file_plot[file_plot$direction_GEO=='-','GEO'] <- 0-file_plot[file_plot$direction_GEO=='-','p_GEO']
file_plot[file_plot$direction_GEO=='+','GEO'] <- file_plot[file_plot$direction_GEO=='+','p_GEO']
file_plot[is.na(file_plot)] <- 0



file_plot=file_plot[,c('celltype','methods','TCGA','GEO')]
tt_plot <- melt(file_plot)
tt_plot[abs(tt_plot$value)<0.1 & tt_plot$value<0 ,'plot_value'] <- -1
tt_plot[abs(tt_plot$value)<0.1 & tt_plot$value>0 ,'plot_value'] <- 1
tt_plot[abs(tt_plot$value)>=0.1 & tt_plot$value<0 ,'plot_value'] <- -0.5
tt_plot[abs(tt_plot$value)>=0.1 & tt_plot$value>0 ,'plot_value'] <- 0.5
tt_plot[is.na(tt_plot)] <- 0
tt_plot$plot_value2 <- tt_plot$plot_value
tt_plot$plot_value <- as.factor(tt_plot$plot_value)


tt_plot$dataset <- tt_plot$methods
tt_plot[tt_plot$methods=='Bayesprism' & tt_plot$variable=='TCGA','dataset'] <- 'Bayesprism2'
tt_plot[tt_plot$methods=='bisque' & tt_plot$variable=='TCGA','dataset'] <- 'bisque2'
tt_plot[tt_plot$methods=='music' & tt_plot$variable=='TCGA','dataset'] <- 'music2'
tt_plot[tt_plot$methods=='DWLS' & tt_plot$variable=='TCGA','dataset'] <- 'DWLS2'
tt_plot[tt_plot$methods=='ReCIDE' & tt_plot$variable=='TCGA','dataset'] <- 'ReCIDE2'
tt_plot[tt_plot$methods=='CIBERSORT' & tt_plot$variable=='TCGA','dataset'] <- 'CIBERSORT2'


tt_plot$dataset <- factor(tt_plot$dataset,c(c('ReCIDE','ReCIDE2', 'music', 'music2','bisque','bisque2', 'Bayesprism','Bayesprism2','DWLS','DWLS2','CIBERSORT','CIBERSORT2')))
tt_plot <- tt_plot[tt_plot$methods=='music',] ##method to select
order_df <- tt_plot[tt_plot$variable=='GEO',]
df <- order_df %>%
  mutate(plot_value = factor(plot_value, levels = c(1, 0.5, -0.5, -1)))


df_sorted <- df %>%
  arrange(plot_value)
celltype_order <- df_sorted[,'celltype']
tt_plot$celltype <- factor(tt_plot$celltype,levels = celltype_order)

main <-  ggplot(tt_plot, aes(x = celltype, y = dataset)) +
  geom_tile(aes(fill = plot_value), color = "#F5F5F5", size = 0.1) +  
  coord_equal() +  
  scale_fill_manual(values = c("-1"="#C1DCFF", "1"="#FFD8E8", "0.5" = "white","-0.5" = "white","0" = "#91C79A")) +  
  geom_text(aes(label = ifelse(abs(plot_value2) == 0.5, ifelse(plot_value2 > 0, "+", "-"), "")), color = 'black', size = 4) +  
  theme(
    legend.position = 'none',
    axis.text.y = element_text(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)
  ) +
  guides(
    fill = guide_colorbar(
      barwidth = 1,        
      barheight = 3,        
      title.position = "top", 
      title.hjust = 0.5    
    )
  )+annotate("rect", xmin = 0.5, xmax = length(unique(tt_plot$celltype)) + 0.5,
             ymin = 0.5, ymax = length(unique(tt_plot$methods))*2 + 0.5,
             color = "black", fill = NA, size = 0.7)


datast_annotation <- tt_plot[,c('dataset','variable')] %>% as.data.frame()
datast_annotation[,'type'] <- 'dataset'
Dataset <-  ggplot(datast_annotation,aes(x=type,y=dataset,fill=variable))+
  geom_tile() +
  coord_equal()+
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  scale_fill_manual(values = c("GEO" ="#fce38a", "TCGA" = "#95e1d3"))+
  theme(legend.position = 'none',panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y= element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5))+
  labs(fill = "Dataset")



method_annotation <- tt_plot[,c('methods','dataset')] %>% as.data.frame()
method_annotation[,'type'] <- 'algorithm'

Method <-  ggplot(method_annotation,aes(x=type,y=methods,fill=dataset))+
  geom_tile() +
  coord_equal()+
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(legend.position = 'left',panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y= element_blank(),axis.text.x = element_blank(),plot.title = element_text(hjust = 0.5))+
  labs(fill = "Method")+scale_fill_manual(values = c("CIBERSORT" ='#ffa500', "DWLS" = "#a08ac1", "Bayesprism" = '#b5d98f', 
                                                     "ReCIDE" = '#e79bbd', "music" ='#aad8f6', "bisque" = '#F58840'))

Method+Dataset+main
main

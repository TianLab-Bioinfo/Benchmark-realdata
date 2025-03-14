setwd('example_dir')  
source('./code/plot/function_figure3b.R')
####################bulk2
ER_files = sort(list.files('example_dir/COAD/GEO/normal/'))
TNBC_files = sort(list.files('example_dir/COAD/GEO/tumor/'))

ER_dir = 'example_dir/COAD/GEO/normal/'
TNBC_dir = 'example_dir/COAD/GEO/tumor/'

key_trend_dir <- "/home/syq/DEC_benmark/TCGA/processed/benchmark_syq/deconvolution_results_final/scRNA_standard/median_total/COAD.rds"
key_trend <- readRDS(key_trend_dir)
row.names(key_trend) <- key_trend[,1]

df=trend_plot(ER_files,ER_dir,TNBC_dir,key_trend,'GEO-COAD')
df=df[df$orig!='scRNA',]

########################bulk1
ER_files = sort(list.files('example_dir/COAD/TCGA/normal/'))
TNBC_files = sort(list.files('example_dir/COAD/TCGA/tumor/'))

ER_dir = 'example_dir/COAD/TCGA/normal/'
TNBC_dir = 'example_dir/COAD/TCGA/tumor/'

key_trend_dir <- "/home/syq/DEC_benmark/TCGA/processed/benchmark_syq/deconvolution_results_final/scRNA_standard/median_total/COAD.rds"
key_trend <- readRDS(key_trend_dir)
row.names(key_trend) <- key_trend[,1]
df2=trend_plot(ER_files,ER_dir,TNBC_dir,key_trend,'TCGA-COAD')

df[is.na(df)] <- 2
df2[is.na(df2)] <- 2

df[,'dataset'] <- 'Bulk1'
df2[,'dataset'] <- 'Bulk2'
plot_df <- rbind(df,df2) %>% data.frame()


method_order <- c(c('ReCIDE','ReCIDE2', 'music', 'music2','bisque','bisque2', 'bayes','bayes2','DWLS','DWLS2','CIBERSORT','CIBERSORT2'))
plot_df$method <- factor(plot_df$orig ,levels = method_order)

order_df <- plot_df[plot_df$orig=='scRNA',]
plot_df <- plot_df[plot_df$orig=='ReCIDE',] ##method to select
plot_df[plot_df$dataset=='Bulk2','orig']<- 'ReCIDE2'

order_df[(order_df$value<0) & abs(order_df$value)<0.05,'plot_value'] <- -1
order_df[(order_df$value>0) & abs(order_df$value)<0.05,'plot_value'] <- 1
order_df[(order_df$value<0) & abs(order_df$value)>0.05,'plot_value'] <- -0.5
order_df[(order_df$value>0) & abs(order_df$value)>0.05,'plot_value'] <- 0.5
order_df[(order_df$value==0) ,'plot_value'] <- 0.3
order_df[(order_df$value==2) ,'plot_value'] <- 0.8
df <- order_df %>%
  mutate(plot_value = factor(plot_value, levels = c(1,0.5, -0.5, -1,0.3,0.8)))


df_sorted <- df %>%
  arrange(plot_value)
celltype_order <- df_sorted[order(df_sorted$plot_value,df_sorted$value),'variable']
plot_df$variable <- factor(plot_df$variable,levels = celltype_order)


plot_df[(plot_df$value<0) & abs(plot_df$value)<0.05,'plot_value'] <- -1
plot_df[(plot_df$value>0) & abs(plot_df$value)<0.05,'plot_value'] <- 1
plot_df[(plot_df$value<0) & abs(plot_df$value)>0.05,'plot_value'] <- -0.5
plot_df[(plot_df$value>0) & abs(plot_df$value)>0.05,'plot_value'] <- 0.5
plot_df[(plot_df$value==0) ,'plot_value'] <- 0.3
plot_df[(plot_df$value==2) ,'plot_value'] <- 0.8

plot_df$plot_value2 <- as.numeric(plot_df$plot_value)
plot_df$plot_value <- as.factor(plot_df$plot_value)

Main <- ggplot(plot_df, aes(x = variable, y = orig)) +
  geom_tile(aes(fill = plot_value), color = "#F5F5F5", size = 0.1) +  
  coord_equal() +  
  scale_fill_manual(values = c("-1"="#C1DCFF", "1"="#FFD8E8", "0.5" = "white","-0.5" = "white","0.8" = "white","0.3" = "#91C79A")) + 
  geom_text(aes(label = ifelse(abs(plot_value2) == 0.5, ifelse(plot_value2 > 0, "+", "-"), "")), color = 'black', size = 4) +  
  theme(
    legend.position = 'none',
    axis.text.y = element_blank(),
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
  )+annotate("rect", xmin = 0.5, xmax = length(unique(plot_df$variable)) + 0.5,
             ymin = 0.5, ymax = length(unique(plot_df$method))*2 + 0.5,
             color = "black", fill = NA, size = 0.7)
Main


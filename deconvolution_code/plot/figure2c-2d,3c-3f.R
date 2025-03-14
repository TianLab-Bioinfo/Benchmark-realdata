library(tidyverse)
library(ggplot2)
library(dplyr)
library(cowplot)



melt_final_df <- readRDS("example_dir/melt_final_df_total.rds")##select final files under conservative or permissive criteria

miss_data <- melt_final_df[melt_final_df$variable=='F1',]## select metrics 
miss_data <- miss_data[!miss_data$disease %in% c('TCGA-PRAD','GEO-PRAD'),]
miss_data$value[is.na(miss_data$value)] <- 0
miss_data$value <- as.numeric(miss_data$value)
mean_values <- data.frame()
for (method in unique(miss_data$method)){mean_values[method,'mean_value'] <- median(miss_data[miss_data$method==method,'value'])}
mean_values["music",'color'] <- '#aad8f6'
mean_values["DWLS",'color'] <- "#a08ac1"
mean_values["CIBERSORT",'color'] <- '#ffa500'
mean_values["bisque",'color'] <- '#F58840'
mean_values["bayes",'color'] <- '#b5d98f'
mean_values["ReCIDE",'color'] <- '#e79bbd'

miss_data[,'method']=factor(miss_data[,'method'],levels=row.names(mean_values[order(mean_values$mean_value),]))
miss_data[,'rank'] <- ave(-miss_data$value, miss_data$disease, FUN = function(x) rank(x, ties.method = "min"))
miss_data[miss_data[,'rank']==1,'significant'] <- TRUE
miss_data[miss_data[,'rank']!=1,'significant'] <- FALSE

vec_color <-mean_values[order(mean_values$mean_value),'color']
ggplot(miss_data,aes(x = method,y = value,fill = method)) +
  geom_boxplot(aes(color = method),
               width = .6, size = .7, alpha = .5 ) +
  theme_classic() +
  scale_color_manual(values= vec_color) +
  scale_fill_manual(values= vec_color) +
  theme(
    axis.text.x = element_text(size = 10, face = "plain", angle = -45),
    axis.text.y = element_text(size = 10, face = "plain"),
    # axis.text.x = element_blank(),
    # axis.text.y = element_blank(),
    axis.title = element_text(size = 8, face = "plain"),
    plot.title = element_text(size = 8, face = "plain", hjust = 0.5),
    plot.subtitle = element_text(size = 10, face = "plain", hjust = 0.5),
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    legend.text = element_text(size = 10),
    legend.spacing.x = unit(0.4, "cm"),
    # axis.title = element_text(size = 8)
    legend.position = "bottom"
  )+
  coord_flip()



miss_data$source <- ifelse(grepl("TCGA", miss_data$disease), "TCGA", "GEO")

miss_data[,'method']=factor(miss_data[,'method'],levels=row.names(mean_values[order(mean_values$mean_value),]))##level和boxplot保持一致
miss_data$source <- factor(miss_data$source, levels = c("TCGA","GEO"))
miss_data <- miss_data %>%
  arrange(source)
miss_data$disease <- factor(miss_data$disease, levels = unique(miss_data$disease))

heatmap_plot <- ggplot(miss_data, aes(x = disease, y = method, fill = value)) +
  geom_tile(color = NA) +  
  scale_fill_gradientn(colors =c('white',"#fed486","#cc0000"), values = c(0, 0.5, 1), limits = c(0, 1)) +
  # labs(title = "Heatmap of Algorithm Values by Dataset", fill = "Value") +
  theme_minimal() +
  #geom_text(aes(label = ifelse(significant, "*", "")), color = "black", size = 3) +
  #coord_fixed()+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank()  
  ) +
  coord_fixed(ratio = 1.3)+
  annotate("rect", xmin = 0.5, xmax = length(unique(miss_data$disease)) + 0.5,
           ymin = 0.5, ymax = length(unique(miss_data$method)) + 0.5,
           color = "black", fill = NA, size = 0.7)



legend_plot <- ggplot(miss_data, aes(x = disease, y = 1, fill = source)) +
  geom_tile() +
  #scale_y_discrete(limits = unique(miss_data$disease))+
  scale_fill_manual(values = c("TCGA" = "#95e1d3", "GEO" = "#fce38a" )) +#
  theme_void() +  
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, -0.5, 0), "cm")  # 减少条栏和热图间距
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.14, 0.14))) +
  coord_fixed(ratio = 0.6)  

final_plot <- plot_grid(legend_plot, heatmap_plot, ncol = 1, rel_heights = c(0.1, 1))


heatmap_plot
ggplot(miss_data,aes(x = method,y = value,fill = method)) +
  geom_boxplot(aes(color = method),#这里的fill如果不设就是空心的
               width = .6, size = .7, alpha = .5, outlier.size = 1, ) +
  theme_classic() +
  scale_color_manual(values= vec_color) +
  scale_fill_manual(values= vec_color) +
  theme(
    axis.text.x = element_text(size = 10, face = "plain", angle = -45),
    axis.text.y = element_text(size = 10, face = "plain"),
    # axis.text.x = element_blank(),
    # axis.text.y = element_blank(),
    axis.title = element_text(size = 8, face = "plain"),
    plot.title = element_text(size = 8, face = "plain", hjust = 0.5),
    plot.subtitle = element_text(size = 10, face = "plain", hjust = 0.5),
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    legend.text = element_text(size = 10),
    legend.spacing.x = unit(0.4, "cm"),
    # axis.title = element_text(size = 8)
    legend.position = "bottom"
  )+stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), vjust = -0.5, size = 3)+
  coord_flip()

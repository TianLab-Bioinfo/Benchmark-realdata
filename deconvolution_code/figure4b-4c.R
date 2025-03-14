
setwd('example_dir/')
library(dplyr)



match_nostrict_confidence  <- readRDS("./data/match_strict.rds")
colnames(match_nostrict_confidence)= c('PAAD','LUSC','LUAD','COAD','HNSC_neg')
miss_data <- melt(as.matrix(match_nostrict_confidence))
colnames(miss_data) <- c('Var1','Var2','value')


miss_data$value <- as.numeric(miss_data$value)
mean_values <- data.frame()
for (method in unique(miss_data$Var1     )){mean_values[method,'mean_value'] <- mean(miss_data[miss_data$Var1==method,'value'])}
mean_values["music",'color'] <- '#aad8f6'
mean_values["DWLS",'color'] <- "#a08ac1"
mean_values["CIBERSORT",'color'] <- '#ffa500'
mean_values["bisque",'color'] <- '#F58840'
mean_values["Bayesprism",'color'] <- '#b5d98f'
mean_values["ReCIDE",'color'] <- '#e79bbd'

miss_data[,'Var1']=factor(miss_data[,'Var1'],levels=row.names(mean_values[order(mean_values$mean_value),]))
miss_data[,'rank'] <- ave(-miss_data$value, miss_data$Var2, FUN = function(x) rank(x, ties.method = "min"))
miss_data[miss_data[,'rank']==1,'significant'] <- TRUE
miss_data[miss_data[,'rank']!=1,'significant'] <- FALSE

vec_color <-mean_values[order(mean_values$mean_value),'color']
ggplot(miss_data,aes(x = Var1,y = value,fill = Var1)) +
  geom_bar(aes(color = Var1),#这里的fill如果不设就是空心的
           stat='summary', fun='mean') +
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
  coord_flip()+
  # scale_y_continuous(
  #   limits = c(0, 1)
  # )+
  coord_cartesian(ylim = c(0,0.8))+
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 3)), vjust = -0.5, size = 3)


# 绘制热图
ggplot(miss_data, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = NA) +  # 去除内部单元格边框
  # '#7DB3D5'
  scale_fill_gradientn(colors =c('white',"#fed486","#cc0000"), values = c(0, 0.5, 1), limits = c(0, 1)) +
  # labs(title = "Heatmap of Algorithm Values by Dataset", fill = "Value") +
  theme_minimal() + 
  #geom_text(aes(label = ifelse(significant, "*", "")), color = "black", size = 3) +
  #coord_fixed()+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank()  # 去除背景网格线
  ) +
  coord_fixed(ratio =0.8)+
  # 添加外部大框
  annotate("rect", xmin = 0.5, xmax = length(unique(miss_data$Var2)) + 0.5,
           ymin = 0.5, ymax = length(unique(miss_data$Var1)) + 0.5,
           color = "black", fill = NA, size = 0.7)
  
  
trend_plot <- function(ER_files,ER_dir,TNBC_dir,key_trend,disease){
  results_list = list()
  disease_type <-  str_split(disease,'-')[[1]][2]
  results_list = list()
  
  for (i in 4:length(ER_files)) {
    if(ER_files[i] == "results_CIBERSORT.txt"){
      CIBERSORT_ER = as.data.frame(read_table(paste0(ER_dir,"results_CIBERSORT.txt"), col_names = TRUE))
      row.names(CIBERSORT_ER) = CIBERSORT_ER[,1]
      CIBERSORT_ER = CIBERSORT_ER[,-1]
      CIBERSORT_ER = CIBERSORT_ER[,-which(colnames(CIBERSORT_ER) %in% c('P-value','Correlation','RMSE'))]
      CIBERSORT_ER = as.data.frame(t(CIBERSORT_ER))
      
      CIBERSORT_TNBC = as.data.frame(read_table(paste0(TNBC_dir,"results_CIBERSORT.txt"), col_names = TRUE))
      row.names(CIBERSORT_TNBC) = CIBERSORT_TNBC[,1]
      CIBERSORT_TNBC = CIBERSORT_TNBC[,-1]
      CIBERSORT_TNBC = CIBERSORT_TNBC[,-which(colnames(CIBERSORT_TNBC) %in% c('P-value','Correlation','RMSE'))]
      CIBERSORT_TNBC = as.data.frame(t(CIBERSORT_TNBC))
      
      results_list[['CIBERSORT']] = list(CIBERSORT_ER,CIBERSORT_TNBC)
    }
    
    if(ER_files[i] == "results_DWLS.rds"){
      DWLS_ER = as.data.frame(readRDS(paste0(ER_dir,"results_DWLS.rds")))
      DWLS_TNBC = as.data.frame(readRDS(paste0(TNBC_dir,"results_DWLS.rds")))
      
      results_list[['DWLS']] = list(DWLS_ER,DWLS_TNBC)
      
    }
    
    if(ER_files[i] == "results_ReCIDE.rds"){
      ReCIDE_ER = as.data.frame(readRDS(paste0(ER_dir,"results_ReCIDE.rds"))[2])
      ReCIDE_TNBC = as.data.frame(readRDS(paste0(TNBC_dir,"results_ReCIDE.rds"))[2])
      
      results_list[['ReCIDE']] = list(ReCIDE_ER,ReCIDE_TNBC)
      
    }
    
    if(ER_files[i] == "results_bayes.rds"){
      bayes_ER = as.data.frame(t(readRDS(paste0(ER_dir,"results_bayes.rds"))))
      bayes_TNBC = as.data.frame(t(readRDS(paste0(TNBC_dir,"results_bayes.rds"))))
      
      results_list[['bayes']] = list(bayes_ER,bayes_TNBC)
      
    }
    
    if(ER_files[i] == "results_bisque.rds"){
      bisque_ER = as.data.frame(readRDS(paste0(ER_dir,"results_bisque.rds")))
      bisque_TNBC = as.data.frame(readRDS(paste0(TNBC_dir,"results_bisque.rds")))
      
      results_list[['bisque']] = list(bisque_ER,bisque_TNBC)
      
    }
    
    if(ER_files[i] == "results_music.rds"){
      music_ER = as.data.frame(t(readRDS(paste0(ER_dir,"results_music.rds"))))
      music_TNBC = as.data.frame(t(readRDS(paste0(TNBC_dir,"results_music.rds"))))
      
      results_list[['music']] = list(music_ER,music_TNBC)
      
    }
  }
  
  if (disease_type =='PRAD'){
    produce_Epi <- function(df){
      name=c("Tumor","Epitheial_Basal","Epitheial_Club","Epitheial_Hillock","Epitheial_Luminal")
      df['Epithelium',] <- apply(df[intersect(name,row.names(df)),],2,sum)
      return(df)
    }
    results_list <-  lapply(results_list, function(pair) {
      lapply(pair, produce_Epi)
    })}
  if (disease_type=='PAAD'){
    produce_Epi <- function(df){
      name=c('Acinar','EMT.Duct','HSP.Duct','Metaplastic','Neoplastic','Normal.Duct')
      df['Epithelium',] <- apply(df[intersect(name,row.names(df)),],2,sum)
      return(df)
    }
    results_list <-  lapply(results_list, function(pair) {
      lapply(pair, produce_Epi)
    })
  }
  
  values_to_keep <-row.names(key_trend)
  
  normalize_df <- function(df){
    df_normalized <- df
    for (i in 1:ncol(df)) {
      row_sum <- sum(df[,i ], na.rm = TRUE)
      if (row_sum != 0) {
        df_normalized[,i ] <- df[,i ] / row_sum
      } else {
        df_normalized[,i ] <- 0  
      }
    }
    
    return(df_normalized)
  }
  
  
  
  keep_rows_with_values <- function(df) {
    df=df[row.names(df) %in% values_to_keep, ]
    return(df)
  }
  
  remove_zero_columns <- function(df) {
    
    non_zero_cols <- sapply(df, function(col) any(col != 0, na.rm = TRUE))
    
    
    df=df[, non_zero_cols, drop = FALSE]
    return(df)
  }
  
  
  
  results_list <- lapply(results_list, function(pair) {
    lapply(pair, keep_rows_with_values)
  })
  results_list <- lapply(results_list, function(pair2) {
    lapply(pair2, normalize_df)
  })
  
  
  for (i in 1:length(results_list)){
    df_category1 = as.data.frame(t(results_list[[i]][[1]]))
    df_category2 = as.data.frame(t(results_list[[i]][[2]]))
    
    col_inter = intersect(colnames(df_category1),colnames(df_category2))
    diff1 <- setdiff(colnames(df_category1), col_inter)
    diff2 <- setdiff(colnames(df_category2), col_inter)
    if(length(diff1) > 0){df_category2[,diff1] = 0}
    
    if(length(diff2) > 0){df_category1[,diff2] = 0}
    
    df_category1 = df_category1[,sort(colnames(df_category1))]
    df_category2 = df_category2[,sort(colnames(df_category2))]
    
    p_values <- numeric(length = ncol(df_category1))
    trends <- character(length = ncol(df_category2))
    diff <- numeric(length = ncol(df_category2))
    
    for (j in seq_along(colnames(df_category1))) {
      col_name <- colnames(df_category1)[j]
      df_category1_in = df_category1[, col_name]
      
      
      df_category2_in = df_category2[, col_name]
      
      
      
      pvalue <- wilcox.test(df_category1_in, df_category2_in, alternative = "two.sided")$p.value
      p_values[j] <- pvalue
      diff[j] <- abs(median(df_category1_in)-median(df_category2_in))
      
      
      if (median(df_category1_in) < median(df_category2_in)) {
        trends[j] <- '+'  ##normal smaller than tumor,tumor high
      } else if (median(df_category1_in) > median(df_category2_in)) {
        trends[j] <- '-'
      } else {
        trends[j] <- '0'  # ?????????????????????
      }
    }
    


    re_df <- data.frame(
      Celltype = colnames(df_category1),
      PValue = adjusted_p_values,
      Trend = trends,
      Diff=diff
    )
    
    re_df$orig=names(results_list)[i]
    if(i == 1){
      RE_DF = re_df
    }else{
      RE_DF <- rbind(RE_DF,re_df)}
  }
  
  
  RE_DF[is.na(RE_DF[,2]),2] = 1
  
  RE_DF[RE_DF$Trend=='-','adj_pvalue'] <- 0- RE_DF[RE_DF$Trend=='-','PValue']
  RE_DF[RE_DF$Trend=='+','adj_pvalue'] <- RE_DF[RE_DF$Trend=='+','PValue']
  RE_DF[RE_DF$Trend==0 & RE_DF$PValue<0.05,'adj_pvalue'] <- 0
  
  
  df=RE_DF[,c('Celltype','orig','adj_pvalue')]
  
  df_wide <- df %>%
    pivot_wider(names_from = Celltype, values_from = adj_pvalue) %>% data.frame()
  melt_df_wide <- melt(df_wide)
  
  
  aa =1
  for (methods in c("CIBERSORT","DWLS","ReCIDE","bayes","bisque","music")) {
    RE_DF_in = RE_DF[RE_DF[,'orig'] == methods,]
    for (mm in 1:nrow(RE_DF_in)) {
      RE_DF_in[mm, 'category'] <- case_when(
        is.na(RE_DF_in[mm, 'Trend']) | is.na(RE_DF_in[mm, 'PValue']) ~ NA_character_,
        RE_DF_in[mm, 'Trend'] == '-' & RE_DF_in[mm, 'PValue'] < 0.05 ~ 'down',
        RE_DF_in[mm, 'Trend'] == '+' & RE_DF_in[mm, 'PValue'] < 0.05 ~ 'up',
        RE_DF_in[mm, 'Trend'] == 0  & RE_DF_in[mm, 'PValue'] < 0.05 ~ 'sig',
        # RE_DF_in[mm, 'Trend'] == '+' & RE_DF_in[mm, 'PValue'] >= 0.05 ~ 'up_no',
        # RE_DF_in[mm, 'Trend'] == '-' & RE_DF_in[mm, 'PValue'] >= 0.05 ~ 'down_no',
        TRUE ~ 'no'
      )
    }
    row.names(RE_DF_in) = RE_DF_in[,1]
    
    RE_DF_in = RE_DF_in[,c('Celltype','orig','category')]
    if(aa == 1){
      p_df = RE_DF_in
      aa = 100
    }else{
      p_df = rbind(p_df,RE_DF_in)
    }
  }
  
  
  plot_df = matrix(nrow = length(unique((p_df[,1]))),
                   ncol = length(unique((p_df[,2]))))
  plot_df = as.data.frame(plot_df)
  row.names(plot_df) = unique((p_df[,1]))
  colnames(plot_df) = unique((p_df[,2]))
  
  
  for (kk in 1:nrow(p_df)) {
    if (!is.na(p_df[kk,3])){
      if(p_df[kk,3] == 'up'){
        plot_df[p_df[kk,1],p_df[kk,2]] = 1}
      if(p_df[kk,3] == 'down'){
        plot_df[p_df[kk,1],p_df[kk,2]] = -1}
      if(p_df[kk,3] == 'no'){
        plot_df[p_df[kk,1],p_df[kk,2]] = 0}
      if(p_df[kk,3] == 'sig'){
        plot_df[p_df[kk,1],p_df[kk,2]] = 0.3}
      if(p_df[kk,3] == 'up_no'){
        plot_df[p_df[kk,1],p_df[kk,2]] = 0.5}
      if(p_df[kk,3] == 'down_no'){
        plot_df[p_df[kk,1],p_df[kk,2]] = -0.5}
    }
    else{p_df[kk,3]  <-  NA}
  }
  
  
  plot_df[] <- lapply(plot_df, as.numeric)
  key_trend[(key_trend$PValue <0.05 )&(key_trend$Trend=='-'),'trend'] <- -1
  key_trend[(key_trend$PValue <0.05 )& (key_trend$Trend=='+'),'trend'] <- 1
  key_trend[(key_trend$PValue > 0.05) & (key_trend$Trend=='+') ,'trend'] <- 0.5
  key_trend[(key_trend$PValue > 0.05) & (key_trend$Trend=='-') ,'trend'] <- -0.5
  
  plot_df <- cbind(plot_df,key_trend[row.names(plot_df),'trend'])
  colnames(plot_df)[ncol(plot_df)]  <- 'Ground_True'
  return(plot_df)
  
}

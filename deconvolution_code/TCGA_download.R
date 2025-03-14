library(dplyr)
library(SingleCellExperiment)


# 查看TCGA中33种癌症的简称
library(TCGAbiolinks)

projects <- TCGAbiolinks::getGDCprojects()$project_id ##获取癌症名字
projects <- projects[grepl('^TCGA', projects, perl=TRUE)]

setwd('./dir_for_TCGA_data/')
for (i in projects){
  clinical_data <- GDCquery_clinic(project = i, type = "clinical")
  dir.create(paste('./',i,sep=''),recursive = T)
  saveRDS(clinical_data,paste('./',i,'/clinical.rds',sep=''))
}


for (j  in projects){
  ###表达数据下载
  query <- GDCquery(
    project = j,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification", 
    workflow.type = "STAR - Counts"
  )
  GDCdownload(query = query)
  expData<- GDCprepare(query = query,
                       save = TRUE
                       
  )
  
  
  df=assay(expData) %>% as.data.frame()
  
  library(fastSave)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  df$ensemble_id=unlist(str_split(row.names(df),"[.]",simplify=T))[,1]
  eg <- bitr(df$ensemble_id,fromType ='ENSEMBL',
             toType = 'SYMBOL',
             OrgDb='org.Hs.eg.db'
  )
  
  dup_name=eg$ENSEMBL[duplicated(eg$ENSEMBL)]
  eg_filter <- eg %>%
    filter(!eg$ENSEMBL %in% dup_name)
  
  for (i in df$ensemble_id){
    if (i %in% eg_filter$ENSEMBL){
      df[df$ensemble_id==i,'Gene'] <- eg_filter[eg_filter$ENSEMBL==i,'SYMBOL']
    }
  }
  
  
  dup_gene=df$Gene[duplicated(df$Gene)]
  df_exp_fil <- df %>%
    filter(!df$Gene %in% dup_gene)
  
  row.names(df_exp_fil) <- df_exp_fil[,ncol(df_exp_fil)]
  df_exp_fil <- df_exp_fil[,-c(1,ncol(df_exp_fil),ncol(df_exp_fil)-1)]
  dir.create(paste('./dir_for_TCGA_data',j,'/Bulk/',sep=''),recursive = T)
  saveRDS.lbzip2(df_exp_fil,file=paste('./dir_for_TCGA_data',j,'/Bulk/total.rdsFS',sep=''),n.cores=30)
  
  
}


counts=seurat_ref@assays[["RNA"]]@counts

eg <- bitr(row.names(counts),fromType ='ENSEMBL',
           toType = 'SYMBOL',
           OrgDb='org.Hs.eg.db'
)

dup_name=eg$ENSEMBL[duplicated(eg$ENSEMBL)]
eg_filter <- eg %>%
  filter(!eg$ENSEMBL %in% dup_name)
row.names(eg_filter) <- eg_filter[,1]
counts_df=merge(counts,eg_filter,by='row.names')
dup_gene=counts_df$SYMBOL[duplicated(counts_df$SYMBOL)]
df_exp_fil <- counts_df %>%
  filter(!counts_df$SYMBOL %in% dup_gene)
row.names(df_exp_fil) <- df_exp_fil[,ncol(df_exp_fil)]
df_exp_fil=df_exp_fil[,-c(1,ncol(df_exp_fil)-1,ncol(df_exp_fil))]


#####download survival information
for (project in projects){
  df_dir=paste('https://gdc-hub.s3.us-east-1.amazonaws.com/download/',project,'.survival.tsv',sep='')
  output_dir=paste('/home/syq/DEC_benmark/TCGA/processed/survival/',project,'_survial.tsv',sep='')
  download.file(df_dir,output_dir)
}

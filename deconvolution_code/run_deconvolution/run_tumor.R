library(Seurat)
library(fastSave)

results_path <- "~/DEC_benmark/TCGA/processed/benchmark_syq/new/COAD/GEO/tumor/"
if (dir.exists(results_path)) {
  cat("Folder exist:", results_path, "\n")
} else {
  cat("Folder not exist")
  dir.create(results_path)
  cat("Folder has been created:", results_path, "\n")
}
source("/home/lmh/ReCIDE_git/ReCIDE/R/ReCIDE_all_function_beta.R")



seurat_ER = readRDS.lbzip2('~/DEC_benmark/TCGA/processed/COAD/scRNA/COAD_tumor.rdsFS',n.cores = 100)
seurat_ER@meta.data[["clMidwayPr"]] <- str_replace_all(seurat_ER@meta.data[["clMidwayPr"]], " ", "_")
seurat_ER@meta.data[["clMidwayPr"]] <- str_replace_all(seurat_ER@meta.data[["clMidwayPr"]], "\\/", "_")
seurat_ER@meta.data[["clMidwayPr"]] <- str_replace_all(seurat_ER@meta.data[["clMidwayPr"]], "\\.", "_")



ER_bulk <- readRDS("/home/syq/DEC_benmark/TCGA/processed/COAD/Bulk/GSE39582/exp_tumor.rds")
ER_bulk[] <- lapply(ER_bulk, function(x) as.integer(as.character(x)))


celltype1 = "clMidwayPr"
subject1 = "PatientTypeID"

dir_ref_recide = paste0(results_path,'/ref_ReCIDE.rds')
dir_results_recide = paste0(results_path,'/results_ReCIDE.rds')
dir_results_bisque = paste0(results_path,'/results_bisque.rds')
dir_results_music = paste0(results_path,'/results_music.rds')
dir_results_bayes = paste0(results_path,'/results_bayes.rds')
dir_ref_CIBERSORT = paste0(results_path,'/ref_CIBERSORT.rds')
dir_results_CIBERSORT = paste0(results_path,'/results_CIBERSORT.txt')
dir_ref_DWLS = paste0(results_path,'/ref_DWLS.rds')
dir_results_DWLS = paste0(results_path,'/results_DWLS.rds')

######run_ReCIDE

source('/home/lmh/ReCIDE/benchmark_syq/main_run_all.R')
func_run_ReCIDE(SC_ref = seurat_ER,
                EXP_df = ER_bulk,
                celltype = celltype1,
                subject = subject1,
                dir_ref = dir_ref_recide,
                dir_results = dir_results_recide)


######run_bisque
source('/home/lmh/ReCIDE/benchmark_syq/main_run_all.R')

func_run_bisque(ref_seurat = seurat_ER,
                bulkdata = ER_bulk,
                celltype = celltype1,
                subject = subject1,
                dir_results = dir_results_bisque)
gc()

######run_music
source('/home/lmh/ReCIDE/benchmark_syq/main_run_all.R')

func_run_music(ref_seurat = seurat_ER,
               EXP_df = ER_bulk,
               celltype = celltype1,
               subject = subject1,
               dir_results = dir_results_music)
gc()




######run_CIBERSORT
source('/home/lmh/ReCIDE/benchmark_syq/main_run_all.R')

func_run_CIBERSORT(combined.data = seurat_ER,
                   bulk.mtx = ER_bulk,
                   celltype = celltype1,
                   subject = subject1,
                   dir_ref = dir_ref_CIBERSORT,
                   dir_results = dir_results_CIBERSORT)
gc()


######run_DWLS
source('/home/lmh/ReCIDE/benchmark_syq/main_run_all.R')

func_run_DWLS(scdata_test = seurat_ER,
              bulkdata = ER_bulk,
              celltype = celltype1,
              dir_DWLS_ref = dir_ref_DWLS,
              dir_DWLS_results = dir_results_DWLS)
gc()


######run_bayes
source('/home/lmh/ReCIDE/benchmark_syq/main_run_all.R')

func_run_bayes(ref_seurat = seurat_ER,
               bulk.mtx = ER_bulk,
               celltype = celltype1,
               dir_results = dir_results_bayes,
               key1 = 'Epi')

# rm(list = ls())
gc()

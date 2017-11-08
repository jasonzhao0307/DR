require(CePa)
require(tidyverse)
require(stringr)
#download raw data from GDAC
Download_GDAC <- function(tumor.name, path.to.folder){
  #STAD no rnaseq2
  tumors.pool = c('LIHC', 'PAAD', 'THCA', 'UCS', 'UCEC','OV','GBM','LGG', 'LIHC', 'BRCA','LUAD','LUSC','PRAD', 'COAD', 'LAML','SKCM','BLCA','STAD','KIRC')
  if (!tumor.name %in% tumors.pool){
    stop(paste0("Please use the name from the list: ", paste(tumors.pool, collapse = ",")))
  }
  gdac.list <- list()
  for (tumor in tumor.name){
    url = paste('http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/', tumor, '/20160128/gdac.broadinstitute.org_', tumor, '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz', sep="")
    if (!file.exists(path.to.folder)){
      dir.create(path.to.folder)
    }
    if (!file.exists(paste0(path.to.folder, '/gdac.broadinstitute.org_', tumor, '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0/', tumor, '.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt'))){
      system(paste0('curl ', '-o ', path.to.folder, '/', tumor,  '.gdac.tar.gz ', url))
      system(paste0('tar -zxvf ', path.to.folder, '/', tumor, '.gdac.tar.gz', ' -C ', path.to.folder))
    }

    df.tmp <- read.table(paste0(path.to.folder, '/gdac.broadinstitute.org_', tumor, '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0/', tumor, '.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt'), stringsAsFactors = F, sep = "\t", header = T, row.names = 1)
    colnames(df.tmp) <- gsub("\\.", "-", colnames(df.tmp))
    df.tmp <- df.tmp[-1,seq(1,ncol(df.tmp),3)]
    df.tmp.numeric <- apply(df.tmp, 2, as.numeric)
    rownames(df.tmp.numeric) <- rownames(df.tmp)
    gdac.list[[paste0(tumor)]] <- df.tmp.numeric
  }
  if (length(gdac.list) == 1){
    return(gdac.list[[1]])
  } else{
  return(gdac.list)
  }
}

#download raw data from GEO
#affy
Download_GEO <- function(geo.accs, path.to.folder){
  geo.list <- list()
  for (geo_acc in geo_accs){
    url = paste('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE2nnn/', geo_acc, '/matrix/', geo_acc, '_series_matrix.txt.gz', sep="")
    if (!file.exists(path.to.folder)){
      dir.create(path.to.folder)
    }

    if (!file.exists(paste0(path.to.folder, "/", geo_acc, '.txt'))){
      system(paste('curl ', url, ' -o ', path.to.folder, '/', geo_acc, '.txt.gz', sep = ""))
      system(paste('gunzip ', path.to.folder, '/', geo_acc, '.txt.gz', sep = ""))
    }
    df.tmp <- read.table(paste0(path.to.folder, "/", geo_acc, '.txt'), stringsAsFactors = F, sep = "\t", header = T, row.names = 1, comment.char = "!")
    geo.list[[paste0(geo_acc)]] <- df.tmp
  }
  if (length(geo.list) == 1){
    return(geo.list[[1]])
  } else{
  return(geo.list)
  }
}


#download raw data from CCLE RPKM
Download_CCLE_RPKM <- function(path.to.folder){
  url = 'https://data.broadinstitute.org/ccle/CCLE_RNAseq_081117.rpkm.gct'

  if (!file.exists(path.to.folder)){
      dir.create(path.to.folder)
  }

  if (!file.exists(paste0(path.to.folder, "/ccle_rpkm.gct"))){
    system(paste('curl ', url, ' -o ', path.to.folder, '/ccle_rpkm.gct', sep = ""))
  }
  df.tmp <- read.gct(paste0(path.to.folder, "/ccle_rpkm.gct"))
  print("read in finished!")
  rowname.vec <- rownames(df.tmp)
  for (i in 1:length(rowname.vec)){
   rowname.vec[i] <- str_match(rowname.vec[i], "(.+)\\.")[2]
  }
  rownames(df.tmp) <- rowname.vec
  return(df.tmp)
}



#download raw data from CCLE metadata
Download_CCLE_Metadata <- function(path.to.folder){
  url = 'https://data.broadinstitute.org/ccle_legacy_data/cell_line_annotations/CCLE_sample_info_file_2012-10-18.txt'
  if (!file.exists(path.to.folder)){
      dir.create(path.to.folder)
  }
  if (!file.exists(paste0(path.to.folder, "/ccle_annotation.txt"))){
    system(paste('curl ', url, ' -o ', path.to.folder, '/ccle_annotation.txt', sep = ""))
  }

  df.tmp <- read.table(paste0(path.to.folder, "/ccle_annotation.txt"), stringsAsFactors = F, sep = "\t", header = T, row.names = 1)
  return(df.tmp)
}

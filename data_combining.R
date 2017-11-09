source("Get_Gene_Mapping_Dict.R")
require(tidyverse)

Map_Data_To_Unifying_Gene_Name <- function(df, database.name){

  # get the dict
  database.name.pool <- c("EXPO", "CCLE", "GDAC")
  gene.name <- ""
  if (database.name == "EXPO"){
    gene.dict <- Get_EXPO_Name_Dict(row.names(df))
    gene.name <- "affy_hg_u133_plus_2"
  } else if (database.name == "CCLE"){
    gene.dict <- Get_CCLE_Name_Dict(row.names(df))
    gene.name <- "ensembl_gene_id"
  } else if (database.name == "GDAC") {
    gene.dict <- Get_GDAC_Name_Dict(row.names(df))
    gene.name <- "entrezgene"
  } else{
    stop(paste0("Please use the name from the list: ", paste(database.name.pool, collapse = ",")))
  }

  gene.dict <- data.frame(apply(gene.dict, 2, as.character), stringsAsFactors = FALSE)

  # index of the original gene id
  index.original.gene <- which(colnames(gene.dict) == gene.name)
  # remove the rows with gene symbol empty
  gene.dict <- gene.dict[which(gene.dict$hgnc_symbol != "" & gene.dict$hgnc_symbol != "?"),]
  #print(head(gene.dict))

  # update the df with only the rows having corresponding gene symbol
  df <- df[which(rownames(df) %in% gene.dict[,index.original.gene]),]

  # get gene symbol names
  gene.symbol.name <- gene.dict$hgnc_symbol[match(rownames(df), gene.dict[,index.original.gene])]



  if (sum(duplicated(gene.symbol.name)) == 0){
    print("No need for aggregation")
    rownames(df) <- gene.symbol.name
    return(df)
  } else{
    df.add.one.col <- cbind.data.frame(df, gene.symbol.name)
    colnames(df.add.one.col)[ncol(df.add.one.col)] <- "GeneID"
    df.add.one.col <- as.tibble(df.add.one.col)
    tmp <- df.add.one.col %>% group_by(GeneID) %>% summarise_all(mean)
    tmp <- as.data.frame(tmp)
    rownames(tmp) <- as.character(tmp$GeneID)
    df <- tmp[,-ncol(tmp)]
    return(df)
  }

}

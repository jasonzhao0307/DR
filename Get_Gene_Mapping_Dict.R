require(biomaRt)
require(stringr)





Get_CCLE_Name_Dict <- function(gene.name){

  for (i in 1:length(gene.name)){
   gene.name[i] <- str_match(gene.name[i], "(.+)\\.")[2]
  }
  dict.output <- Get_Gene_Mapping_Dict("ensembl_gene_id", gene.name)
  return(dict.output)

}


Get_GDAC_Name_Dict <- function(gene.name){
  entrez.vec <- c()
  gene.symbol.vec <- c()
  for (i in 1:length(gene.name)){
    entrez.vec[i] <- str_match(gene.name[i], "\\|([0-9]+)")[2]
    gene.symbol.vec[i] <- str_match(gene.name[i], "(.+)\\|[0-9]+")[2]
  }
  df.output <- data.frame(entrezgene = entrez.vec, hgnc_symbol = gene.symbol.vec)
  return(df.output)
}


Get_EXPO_Name_Dict <- function(gene.name){
  dict.output <- Get_Gene_Mapping_Dict("affy_hg_u133_plus_2", gene.name)
  return(dict.output)
}


# get dictionary for mapping
Get_Gene_Mapping_Dict <- function(filter, gene.name){
  filter.pool <- c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "affy_hg_u133_plus_2")
  if (!filter %in% filter.pool){
    stop(paste0("Please use the name from the list: ", paste(filter.pool, collapse = ",")))
  }
  ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  name.dict <- getBM(attributes = c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "affy_hg_u133_plus_2"), filters = filter, values = gene.name, mart=ensembl)
  #name.dict <- apply(name.dict, 2, as.character)
  return(name.dict)
}



# get the updated expression matrix in gene symbol

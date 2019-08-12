suppressMessages(library(dplyr))
suppressMessages(library(foreach))
suppressMessages(library(doSNOW))
library(xml2)
library(openxlsx)
suppressMessages(library(doParallel))

#===================== FUNCTION =====================

get_info <- function(proid){
  library(xml2)
  library(dplyr)
  name <- paste0("http://www.uniprot.org/uniprot/", proid, ".xml")
  result <- try(read_xml(name), silent = T)
  
  error <- 1
  if (!class(result)[1] == "xml_document"){
    while(any(class(result) == "try-error"& error == 1)){
      cat(paste0("The error protein ID is ", proid, ". If this error happens again, please check your network or protein ID !!!\n"))
      if (grepl("Failed to parse text", result[1])){
        error <- 0
      }else{
        result <- try(read_xml(name), silent = T)
      }
    }
  }
  
  if (any(class(result) == "try-error")){
    out <- c(ProteinID = proid, ProteinName = "-", BP = "-", MF = "-", CC = "-", GO_Term = "-")
  }else{
    # goid <- xpathSApply(result, "//g:dbReference[@type='GO']", xmlGetAttr, "id", namespaces = "g")
    # tmp <- xpathSApply(result, "//j:property [@type='term']", xmlGetAttr, "value", namespaces = "j")
    # category <- unlist(lapply(tmp, function(x){strsplit(x, ":")[[1]][1]}))
    # goterm <- unlist(lapply(tmp, function(x){strsplit(x, ":")[[1]][2]}))
    # proname <- xpathSApply(result, "//j:fullName", fun = xmlValue, namespaces = "j")[1]
    
    goid <- xml_find_all(result, ".//d1:dbReference[@type='GO']", ns = xml_ns(result)) %>% xml_attr("id")
    tmp <- xml_find_all(result, ".//d1:property[@type='term']", ns = xml_ns(result)) %>% xml_attr("value")
    category <- unlist(lapply(tmp, function(x){strsplit(x, ":")[[1]][1]}))
    goterm <- unlist(lapply(tmp, function(x){strsplit(x, ":")[[1]][2]}))
    proname <- xml_find_all(result, ".//d1:fullName", ns = xml_ns(result)) %>% xml_text() %>% head(1)
    geneid <- xml_find_all(result, ".//d1:dbReference[@type='GeneID']", ns = xml_ns(result)) %>% xml_attr("id") %>% head(1)
    genename <- xml_find_all(result, ".//d1:gene/d1:name", ns = xml_ns(result)) %>% xml_text() %>% head(1)
    
    proname <- ifelse(length(proname) == 0, NA, proname)
    geneid <- ifelse(length(geneid) == 0, NA, geneid)
    genename <- ifelse(length(genename) == 0, NA, genename)
    
    if ("P" %in% category){
      bp <- paste0(paste(goterm[which(category == "P")], paste0("[", goid[which(category == "P")], "]")), collapse = ";")
    }else{
      bp <- "-"
    }
    
    if ("F" %in% category){
      mf <- paste0(paste(goterm[which(category == "F")], paste0("[", goid[which(category == "F")], "]")), collapse = ";")
    }else{
      mf <- "-"
    }
    
    if ("C" %in% category){
      cc <- paste0(paste(goterm[which(category == "C")], paste0("[", goid[which(category == "C")], "]")), collapse = ";")
    }else{
      cc <- "-"
    }
    
    all_go <- paste(c(bp[!is.na(bp)], mf[!is.na(mf)], cc[!is.na(cc)]), collapse = ";")
    
    if (length(goid) == 0){
      out <- c(ProteinID = proid, ProteinName = proname, gene_id = geneid, gene_name = genename, BP = "-", MF = "-", CC = "-", GO_Term = "-")
    }else{
      out <- c(ProteinID = proid, ProteinName = proname, gene_id = geneid, gene_name = genename, BP = bp, MF = mf, CC = cc, GO_Term = all_go)
    }
  }
  return(out)
}

uniprot_id <- read.table(file = "../Desktop/protein_list.txt", header = F, stringsAsFactors = F) %>% as.matrix() %>% as.character()

cl.cores <- detectCores()
print(cl.cores)
core <- 4
cl <- makeCluster(core)
registerDoSNOW(cl)
pb <- txtProgressBar(min = 0, max = length(uniprot_id), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

system.time({
  result_df <- foreach(i = 1:length(uniprot_id), .combine = rbind, .options.snow = opts) %dopar% {
    df <- get_info(proid = uniprot_id[i])
    
    Sys.sleep(0.2)
    return(df)
  }
  close(pb)
  stopCluster(cl)
})

write.table(result_df, file = "Uniprot_API_result.xls", sep = "\t", col.names = T, rrow.names = F, quote = F)


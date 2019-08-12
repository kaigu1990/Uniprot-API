library(httr)
library(jsonlite)

#===================== FUNCTION =====================

get_url_info <- function(url){
  d <- GET(url, accept("application/json"))
  mydata <- toJSON(content(d))
  mydata <- fromJSON(mydata)
  
  res <- matrix(nrow = nrow(mydata), ncol = 8)
  for(i in 1:nrow(mydata)){
    proteinid <- mydata$accession[[i]]
    proteinname <- ifelse(!is.null(mydata$protein$recommendedName$fullName$value[[i]]), mydata$protein$recommendedName$fullName$value[[i]],
                          ifelse(!is.null(mydata$protein$submittedName[[i]]$fullName$value[[1]]), mydata$protein$submittedName[[i]]$fullName$value[[1]],
                                 ifelse(!is.null(mydata$protein$alternativeName[[i]]$fullName$value[[1]]), mydata$protein$alternativeName[[i]]$fullName$value[[1]],
                                        NA)))
    genename <- ifelse(!is.null(mydata$gene[[i]]$name$value[[1]]), mydata$gene[[i]]$name$value[[1]],
                       ifelse(!is.null(mydata$gene[[i]]$orfNames[[1]]$value[[1]]), mydata$gene[[i]]$orfNames[[1]]$value[[1]],
                              NA))
    temp <- data.frame(type = unlist(mydata$dbReferences[[i]]$type), id = unlist(mydata$dbReferences[[i]]$id), stringsAsFactors = F)
    geneid <- temp$id[which(temp$type == "GeneID")][1]
    
    go_term <- unlist(mydata$dbReferences[[i]]$properties$term)
    go_id <- temp$id[which(temp$type == "GO")]
    
    go_P_index <- grepl(pattern = "^P:", go_term)
    if (any(go_P_index)){
      go_P_Term <- unlist(lapply(go_term[go_P_index], function(x){
        strsplit(x, "P:")[[1]][2]
      }))
      go_P <- paste0(go_P_Term, " ", "[", go_id[go_P_index], "]")
    }else{
      go_P <- NA
    }
    
    go_F_index <- grepl(pattern = "^F:", go_term)
    if (any(go_F_index)){
      go_F_Term <- unlist(lapply(go_term[go_F_index], function(x){
        strsplit(x, "F:")[[1]][2]
      }))
      go_F <- paste0(go_F_Term, " ", "[", go_id[go_F_index], "]")
    }else{
      go_F <- NA
    }
    
    go_C_index <- grepl(pattern = "^C:", go_term)
    if (any(go_C_index)){
      go_C_Term <- unlist(lapply(go_term[go_C_index], function(x){
        strsplit(x, "C:")[[1]][2]
      }))
      go_C <- paste0(go_C_Term, " ", "[", go_id[go_C_index], "]")
    }else{
      go_C <- NA
    }
    
    go_all <- na.omit(c(go_P, go_F, go_C))
    if (length(go_all) == 0){
      go_all <- NA
    }
    
    res[i,] <- c(proteinid, proteinname, geneid, genename, paste(go_P, collapse = ";"), paste(go_F, collapse = ";"), paste(go_C, collapse = ";"), paste(go_all, collapse = ";"))
  }
  
  res <- data.frame(res, stringsAsFactors = F)
  names(res) <- c("ProteinID", "ProteinName", "gene_id", "gene_name", "BP", "MF", "CC", "GO_Term")
  
  return(res)
}

uniprot_id <- read.table(file = "../Desktop/protein_list.txt", header = F, stringsAsFactors = F) %>% as.matrix() %>% as.character()
# id <- uniprot_id
res <- data.frame(stringsAsFactors = F)
n <- seq(1, length(uniprot_id), by = 100)
for(i in n){
  sub_id <- uniprot_id[i:(i+99)]
  sub_id <- na.omit(sub_id)
  sub_id <- paste(sub_id, collapse = ",")
  url <- paste0("https://www.ebi.ac.uk/proteins/api/proteins?", "accession=", sub_id)
  
  subres <- get_url_info(url = url)
  res <- dplyr::bind_rows(res, subres)
}

write.table(res, file = "Uniprot_API_result.xls", sep = "\t", col.names = T, rrow.names = F, quote = F)


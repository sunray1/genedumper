#this script gets all species on Genbank for a specific taxon
#remotes::install_github("ropensci/taxize")
library(taxize)
library(stringr)
#csvin <- read.csv(file = "Taxon.csv", header=FALSE)
#csvlist <- levels(csvin$V1)

#edit these
Sys.setenv(ENTREZ_KEY="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
csvlist <- "xxxxxxxxxxxxxx"
divisionfilter <- "xxxxxxxxxxxx"
#automatically gets genus and species, don't add
ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "subgenus")

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

ranks <- sapply(ranks, simpleCap)

sp_out <- downstream(csvlist, downto = "Species", db = "ncbi", ambiguous_nodes=TRUE, ambiguous_species=FALSE)
sp_list <- sp_out[[1]]$childtaxa_name
csvlist_split <- str_split_fixed(sp_list, " ", 2)
colnames(csvlist_split) <- c("Genus", "Epithet")
genera <- unique(csvlist_split[,1])

if ('Subgenus' %in% ranks || 'subgenus' %in% ranks){
  class_test <- tax_name(query = sp_list, get = c(ranks,"Genus"), db = "ncbi", rank_filter = "species", division_filter = divisionfilter)
  taxa_merge <- merge(csvlist_split, class_test, by.x="Genus", by.y="Genus")
  final_tax <- cbind(taxa_merge[,5:ncol(class_test)], taxa_merge[,1], taxa_merge[,ncol(class_test)+1], paste(taxa_merge[,1], taxa_merge[,2], sep = " "), stringsAsFactors = FALSE)
  colnames(final_tax) <- c(ranks[1:match("Subgenus", ranks)-1], "Genus", "Subgenus", "Species")
} else {
  class_test <- tax_name(query = genera, get = ranks, db = "ncbi", rank_filter = "genus", division_filter = divisionfilter)
  taxa_merge <- merge(csvlist_split, class_test, by.x="Genus", by.y="query")
  final_tax <- cbind(taxa_merge[,3:ncol(class_test)+1], taxa_merge[,1], paste(taxa_merge[,1], taxa_merge[,2], sep = " "), stringsAsFactors = FALSE)
  colnames(final_tax) <- c(ranks, c('Genus', 'Species'))
}


#this will change if we change the ranks


final_tax[is.na(final_tax)] <- ""

write.csv(final_tax, file = "Taxonomy.csv", row.names = FALSE, quote=FALSE)

#this script gets all species on Genbank for a specific taxon
#remotes::install_github("ropensci/taxize")
library(taxize)
library(stringr)
Sys.setenv(ENTREZ_KEY="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
#csvin <- read.csv(file = "Taxon.csv", header=FALSE)
#csvlist <- levels(csvin$V1)
csvlist <- "xxxxxx"
#automatically gets genus and species, don't add
ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family")


sp_out <- downstream(csvlist, downto = "Species", db = "ncbi", ambiguous=TRUE)
sp_list <- sp_out[[1]]$childtaxa_name
csvlist_split <- str_split_fixed(sp_list, " ", 2)
genera <- unique(csvlist_split[,1])

class_test <- tax_name(query = genera, get = ranks, db = "ncbi")
#can add division_filter = "moths" to reduce interactiveness
colnames(csvlist_split) <- c("Genus", "Epithet")
taxa_merge <- merge(csvlist_split, class_test, by.x="Genus", by.y="query")

#this will change if we change the ranks
final_tax <- cbind(taxa_merge[,3:ncol(class_test)+1], taxa_merge[,1], paste(taxa_merge[,1], taxa_merge[,2], sep = " "))
colnames(final_tax) <- c(ranks, c('Genus', 'Species'))
final_tax[is.na(final_tax)] <- ""

write.csv(final_tax, file = "Taxonomy.csv", row.names = FALSE, quote=FALSE)

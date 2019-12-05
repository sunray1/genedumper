#this script gets all species on Genbank for a specific taxon

library(taxize)
library(stringr)
Sys.setenv(ENTREZ_KEY="ENTER_ENTREZ_KEY_HERE")
csvin <- read.csv(file = "Taxon.csv", header=FALSE)
csvlist <- levels(csvin$V1)
#csvlist <- "Pieridae"


sp_out <- downstream(csvlist, downto = "Species", db = "ncbi")
sp_list <- sp_out[[1]]$childtaxa_name
csvlist_split <- str_split_fixed(sp_list, " ", 2)
genera <- unique(csvlist_split[,1])

class_test2 <- tax_name(query = genera, get = c("Kingdom", "Phylum", "Class", "Order", 'Family'), db = "ncbi")
#can add division_filter = "moths" to reduce interactiveness
colnames(csvlist_split) <- c("Genus", "Epithet")
taxa_merge <- merge(csvlist_split, class_test, by.x="Genus", by.y="query")
taxise_out <- cbind(taxa_merge[,3:ncol(class_test)],class_test[,2])
final_tax <- cbind(taxa_merge[,4:ncol(taxa_merge)],taxa_merge[,1], paste(taxa_merge[,1], taxa_merge[,2], sep = " "))
colnames(final_tax) <- c("Kingdom", "Phylum", "Class", "Order", 'Family', 'Genus', 'Species')
final_tax[is.na(final_tax)] <- ""

write.csv(final_tax, file = "Taxonomy.csv", row.names = FALSE, quote=FALSE)

library(taxize)
library(stringr)
Sys.setenv(ENTREZ_KEY="2ab4a57550b4801d78fa9ac204b373348308")
#csvin <- read.csv(file = "Species_list.csv", header=FALSE)
#csvlist <- levels(csvin$V1)
csvlist <- "Pieridae"

#this takes a while to run
sp_out <- downstream(csvlist, downto = "Species", db = "ncbi")
sp_list <- sp_out[[1]]$childtaxa_name
csvlist_split <- str_split_fixed(sp_list, " ", 2)
genera <- unique(csvlist_split[,1])

class_test2 <- tax_name(query = genera, get = c("Kingdom", "Phylum", "Class", "Order", 'Family', 'Subfamily', 'Tribe'), db = "ncbi")
#can add division_filter = "moths" or division_filter = "butterflies" to reduce interactiveness to above
colnames(csvlist_split) <- c("Genus", "Epithet")
taxa_merge <- merge(csvlist_split, class_test2, by.x="Genus", by.y="query")
#taxise_out <- cbind(taxa_merge[,3:ncol(class_test2)],class_test2[,2])
final_tax <- cbind(taxa_merge[,4:ncol(taxa_merge)],taxa_merge[,1], paste(taxa_merge[,1], taxa_merge[,2], sep = " "))
colnames(final_tax) <- c("Kingdom", "Phylum", "Class", "Order", 'Family', 'Subfamily', 'Tribe', 'Genus', 'Species')
final_tax[is.na(final_tax)] <- ""

write.csv(final_tax, file = "Taxonomy.csv", row.names = FALSE, quote=FALSE)

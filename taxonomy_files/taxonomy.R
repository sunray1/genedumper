library(taxize)
Sys.setenv(ENTREZ_KEY="2ab4a57550b4801d78fa9ac204b373348308")
csvin <- read.csv(file = "genus_list.csv", header=FALSE)
csvlist <- levels(csvin$V1)
#csv_first <- csvlist[1:10]
#csv_second <- csvlist[423:843]
class_test <- tax_name(query = csvlist, get = c("Kingdom", "Phylum", "Class", "Order", 'Family'), db = "ncbi")
#class_test <- classification(csv_first, db="ncbi")
write.csv(class_test[,2:7], file = "Taxonomy.csv", row.names = FALSE, quote=FALSE)

#class_test <- classification(csvlist, db="gbif")
#for (i in class_test){
#  write(str_c(i[1], collapse = ","), file = "taxonomy.csv", append=TRUE)
#}


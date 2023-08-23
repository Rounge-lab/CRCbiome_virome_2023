#!/usr/bin/Rscript


checkV_format <- do.call("rbind", lapply(snakemake@input[["checkV_files"]], function(f) {
  read.delim(f, stringsAsFactors = FALSE, header = FALSE)
}))

new_headers <- c("Marker lineage", "# genomes", "# markers", "# marker sets", "0", "1", "2", "3", "4", "5+", "Strain heterogeneity")

checkM_format <- do.call("rbind", lapply(snakemake@input[["checkV_files"]], function(f) {
  tmp <- read.delim(f, stringsAsFactors = FALSE, header = FALSE)
  
  # tmp$V1 <- ifelse(tmp$V3 %in% "Yes", paste(tmp$V1, "_1", sep = ""), tmp$V1)
  
  tmp$Completeness <- tmp$V10
  # tmp$Contamination <- tmp$V12
  tmp$Contamination <- 0
  tmp$`Bin Id` <- tmp$V1
  
  tmp <- cbind(tmp, sapply(new_headers, function(x) rep(0,nrow(tmp))))
  tmp[ , c("Bin Id", "Marker lineage", "# genomes", "# markers", "# marker sets", "0", "1", "2", "3", "4", "5+", "Completeness", "Contamination", "Strain heterogeneity")]
}))

write.table(checkV_format, file = snakemake@output[["checkV_gathered"]], sep = "\t", quote = FALSE, row.names = FALSE)
write.table(checkM_format, file = snakemake@output[["checkM_format"]], sep = "\t", quote = FALSE, row.names = FALSE)

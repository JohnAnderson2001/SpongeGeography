
taxonomy <- read.csv("reference_with_taxonomy.csv", stringsAsFactors = FALSE)
taxonomy_pruned <- taxonomy[,c("Species", 'class', 'subclass', 'order', 'family')]
constraint_result <- data.frame()
cat("\n", file="cladeconstraints.Rev", append=FALSE)
final_focal_string <- c()
#for (rank in c("class", "subclass", "order")) {
for (rank in c("order")) {

	taxa <- unique(taxonomy_pruned[,rank])
	taxa <- taxa[!is.na(taxa)]
	taxonomy_pruned[is.na(taxonomy_pruned[,rank]), rank] <- ""
	for (taxon in taxa) {
		focal_species <- taxonomy_pruned[taxonomy_pruned[,rank] == taxon, "Species"]
		if(length(focal_species) > 1) {
			focal_species <- gsub(" ", "_", focal_species)
			focal_species <- paste0('"', focal_species, '"')
			constraint_result <- rbind(constraint_result, data.frame(clade=taxon, species=paste0(focal_species, collapse = " ")))
			cat("clade_", rank, "_", taxon, " = clade(", paste0(focal_species, collapse = ", "), ")\n\n", file="cladeconstraints.Rev", append=TRUE, sep="")	
			final_focal_string <- c(final_focal_string, paste0("clade_", rank, "_", taxon))
		}
	
	}
}
cat(paste0("constraints = v(", paste0(final_focal_string, collapse=", "), ")\n"), file="cladeconstraints.Rev", append=TRUE)

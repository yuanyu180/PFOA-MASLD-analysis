library(dplyr)
library(tidyr)
uniprot_df <- read.table("refer.UniProt.txt", sep = "\t", header = TRUE, fill = TRUE, na.strings = c("", "NA"))
chembl_df <- read.table("input.tsv", sep = "\t", header = TRUE)
chembl_expanded_df <- chembl_df %>%
  separate_rows(UniProt.Accessions, sep = "\\|")
merged_df <- inner_join(uniprot_df, chembl_expanded_df, by = c("Entry" = "UniProt.Accessions"))
merged_df <- merged_df %>%
  mutate(gene_name = gsub("_HUMAN", "", Entry.Name))
print("Merged DataFrame with gene_name:")
print(head(merged_df))
write.table(merged_df, "merged_output.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Given vector of cell line names, get cellosaurus accession ids
library(AnnotationGx)
library(data.table)

file <- read.csv("input_data/unmatched_cell_lines_adjusted.csv")
names <- file$modified_name

# Map cell lines to cellosaurus accessions
samples <- mapCell2Accession(names)

# Get existing cell line accession list and see which compare with new accessions to find any unique
existing_accessions <- read.csv("input_data/cellosaurus_ids.csv")
non_overlap_accessions <- setdiff(samples$accession, existing_accessions$cellosaurus_id)

# conver to df to allow for appending to existing cell lines
new_accessions_df <- data.frame(cellosaurus_id = non_overlap_accessions)

print(new_accessions_df)

write.table(
    new_accessions_df,
    file = "input_data/cellosaurus_ids.csv",
    row.names = FALSE,
    col.names = FALSE,
    quote = TRUE,
    append = TRUE
)



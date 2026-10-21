# Given vector of cell line names, get cellosaurus accession ids
library(AnnotationGx)
library(data.table)

file <- read.csv("input_data/cell_lines_ctrp.csv")
names <- file$name

# Map cell lines to cellosaurus accessions
samples <- mapCell2Accession(names)
samples <- samples[!is.na(samples$accession)]

# Get existing cell line accession list and see which compare with new accessions to find any unique
existing_accessions <- read.csv("input_data/cellosaurus_ids.csv")
non_overlap_accessions <- setdiff(samples$accession, existing_accessions$cellosaurus_id)

# conver to df to allow for appending to existing cell lines
new_accessions_df <- data.frame(cellosaurus_id = non_overlap_accessions)

print(new_accessions_df)

# Ensure existing file ends with a newline before appending
csv_path <- "input_data/cellosaurus_ids.csv"
raw <- readChar(csv_path, file.info(csv_path)$size)
if (nchar(raw) > 0 && !grepl("\n$", raw)) {
    cat("\n", file = csv_path, append = TRUE)
}

write.table(
    new_accessions_df,
    file = csv_path,
    row.names = FALSE,
    col.names = FALSE,
    quote = TRUE,
    append = TRUE
)

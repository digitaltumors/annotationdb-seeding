library(data.table)
f <- "output_data/union/mar-27-2026/union_out.csv"
dt <- fread(f, colClasses = "character")
before <- nrow(dt)
dt <- dt[!is.na(inchikey) & nzchar(trimws(inchikey)) & !is.na(smiles) & nzchar(trimws(smiles))]
after <- nrow(dt)
cat(sprintf("Removed %d rows missing InChIKey or SMILES.\n", before - after))
fwrite(dt, f)
cat(sprintf("Remaining rows: %d\n", after))
cat("File size:", file.info(f)$size / 1e6, "MB\n")

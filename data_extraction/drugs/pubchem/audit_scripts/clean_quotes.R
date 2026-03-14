library(data.table)
f <- "output_data/union/complete/union_out.csv"
dt <- fread(f, colClasses = "character")

# Find lengths before cleanup
lens_before <- sapply(dt, function(x) sum(nchar(x), na.rm = TRUE))
cat("Total chars before:", sum(lens_before), "\n")

# Remove all double quotes from all columns
for (col in names(dt)) {
    dt[, (col) := gsub('"', "", get(col), fixed = TRUE)]
    dt[, (col) := trimws(get(col))]
}

# Find lengths after cleanup
lens_after <- sapply(dt, function(x) sum(nchar(x), na.rm = TRUE))
cat("Total chars after:", sum(lens_after), "\n")

# Re-write the file
fwrite(dt, "output_data/union/complete/union_out_clean.csv")
cat("Clean file written to output_data/union/complete/union_out_clean.csv\n")

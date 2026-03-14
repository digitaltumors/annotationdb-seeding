library(data.table)

f <- "output_data/union/complete/union_out.csv"
dt <- fread(f, colClasses = "character")

before <- nrow(dt)
cat(sprintf("Before dedup: %d rows\n", before))

# For duplicate InChIKeys, merge mapped_name by combining all unique non-empty tokens
dt <- dt[,
    {
        if (.N == 1) {
            .SD
        } else {
            # Merge all mapped_names across duplicate rows
            all_names <- unlist(strsplit(paste(mapped_name, collapse = "; "), ";\\s*"))
            all_names <- unique(trimws(all_names[nzchar(trimws(all_names))]))
            merged_name <- if (length(all_names)) paste(all_names, collapse = "; ") else ""
            # Take the first row as base, update mapped_name
            result <- .SD[1]
            result[, mapped_name := merged_name]
            result
        }
    },
    by = inchikey
]

after <- nrow(dt)
cat(sprintf("After dedup : %d rows (removed %d duplicates)\n", after, before - after))

fwrite(dt, f)
cat("Saved.\n")

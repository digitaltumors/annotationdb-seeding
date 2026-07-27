library(AnnotationGx)
library(data.table)
library(jsonlite)

# Import drug list
cell_lines_file <- read.csv("input_data/cellosaurus_ids.csv")
cell_lines <- cell_lines_file$cellosaurus_id


# Map cell line anmes to cellosaurus accession id, remove NA's
# result <- mapCell2Accession(cell_lines, parsed = TRUE)
# result <- result[!is.na(result$accession)]

flatten_col <- function(x) {
    if (!is.list(x)) return(x)
    vapply(
        x,
        function(el) {
            if (is.null(el) || length(el) == 0) {
                NA_character_
            } else if (is.list(el) && !is.null(names(el))) {
                # Named/nested list: keep structure as JSON
                jsonlite::toJSON(el, auto_unbox = TRUE, null = "null")
            } else {
                # Vector/unnamed list: collapse unique values
                paste(unique(as.character(unlist(el))), collapse = " | ")
            }
        },
        character(1)
    )
}


# List wanted cellosaurus fields
fields <- c(
  "id", "sy", "idsy", "ac", "acas", "dr", "ww", "genome-ancestry", "hla",
  "registration", "sequence-variation", "anecdotal", "biotechnology", "breed",
  "caution", "cell-type", "characteristics", "donor-info", "derived-from-site",
  "discontinued", "doubling-time", "from", "group", "karyotype", "knockout",
  "msi", "miscellaneous", "misspelling", "mab-isotype", "mab-target", "omics",
  "part-of", "population", "problematic", "resistance", "senescence",
  "integrated", "transformant", "virology", "cc", "str", "di", "din", "dio",
  "ox", "sx", "ag", "oi", "hi", "ch", "ca", "dt", "dtc", "dtu", "dtv"
)
field_data <- annotateCellAccession(cell_lines, to = fields)
field_data_flat <- as.data.frame(lapply(field_data, flatten_col), stringsAsFactors = FALSE)

# Build a normalized synonyms table including Cellosaurus synonyms and crossReferences.
# Use ONLY the 'accession' field for mapping.
dt <- as.data.table(field_data)

# --- minimal aliasing so old column names exist (no structural change) ---
if (!("accession" %in% names(dt)) && ("ac" %in% names(dt))) {
  dt[, accession := vapply(ac, function(x) if (is.null(x) || length(x) == 0) NA_character_ else as.character(x[[1]]), character(1))]
}
if (!("synonyms" %in% names(dt)) && ("sy" %in% names(dt))) {
  dt[, synonyms := sy]
}
if (!("crossReferences" %in% names(dt)) && ("dr" %in% names(dt))) {
  dt[, crossReferences := dr]
}
# -----------------------------------------------------------------------

# --- de-duplicate by primary accession BEFORE any exports (also keeps flat table in sync) ---
acc_vec <- dt$accession
valid_idx <- !is.na(acc_vec) & nzchar(acc_vec)
dedup_idx <- valid_idx & !duplicated(acc_vec)
dt <- dt[dedup_idx]
field_data_flat <- field_data_flat[dedup_idx, , drop = FALSE]
# --- END ---

# Ensure required columns exist
if (!("accession" %in% names(dt))) {
    stop("Expected 'accession' column not found. Ensure 'ac' was requested in 'fields'.")
}
if (!("synonyms" %in% names(dt))) {
    dt[, synonyms := vector("list", .N)]
}
if (!("crossReferences" %in% names(dt))) {
    dt[, crossReferences := vector("list", .N)]
}

# 1) Regular Cellosaurus synonyms -> source = "cellosaurus"
syn_cellosaurus <- rbindlist(
    lapply(seq_len(nrow(dt)), function(i) {
        syns <- dt$synonyms[[i]]
        if (is.null(syns) || length(syns) == 0) return(NULL)
        data.table(
            accession = dt$accession[i],
            synonym = trimws(as.character(unlist(syns))),
            source = "cellosaurus"
        )
    }),
    use.names = TRUE, fill = TRUE
)

# 2) Cross-references as synonyms, keeping their source (e.g., GEO, DepMap, etc.)
syn_crossrefs <- rbindlist(
    lapply(seq_len(nrow(dt)), function(i) {
        cr <- dt$crossReferences[[i]]
        if (is.null(cr) || length(cr) == 0) return(NULL)
        # Iterate by index to preserve duplicate names if present
        rbindlist(
            lapply(seq_along(cr), function(j) {
                vals <- cr[[j]]
                if (is.null(vals) || length(vals) == 0) return(NULL)
                src <- names(cr)[j]
                data.table(
                    accession = dt$accession[i],
                    synonym = trimws(as.character(unlist(vals))),
                    source = ifelse(is.null(src) || identical(src, ""), NA_character_, src)
                )
            }),
            use.names = TRUE, fill = TRUE
        )
    }),
    use.names = TRUE, fill = TRUE
)

# Combine, clean
syn_dt <- rbindlist(list(syn_cellosaurus, syn_crossrefs), use.names = TRUE, fill = TRUE)
if (is.null(syn_dt) || nrow(syn_dt) == 0L) {
  syn_dt <- data.table(accession = character(0), synonym = character(0), source = character(0))
}
syn_dt <- syn_dt[!is.na(synonym) & nzchar(synonym)]
# normalize source and synonym (trim + collapse spaces) for consistent handling
syn_dt[, source := tolower(trimws(source))]
syn_dt[, synonym := gsub("\\s+", " ", trimws(synonym), perl = TRUE)]

# --- NEW: collapse to ONE row per (accession, synonym) case-insensitively,
#          still preferring source="cellosaurus"
if (nrow(syn_dt)) {
  syn_dt[, norm_syn := tolower(synonym)]
  syn_dt[, pref := fifelse(source == "cellosaurus", 0L, 1L)]
  setorder(syn_dt, accession, norm_syn, pref, source, synonym)
  syn_dt <- syn_dt[!duplicated(syn_dt, by = c("accession", "norm_syn"))]
  syn_dt[, c("pref", "norm_syn") := NULL]
}
# --- END NEW ---


write.csv(field_data_flat, "output_data/cell_line_table.csv", row.names = TRUE)

# match DB FK/PK column name
if ("accession" %in% names(syn_dt)) {
  data.table::setnames(syn_dt, "accession", "cellosaurus_accession")
}
write.csv(as.data.frame(syn_dt), "output_data/cell_line_synonyms.csv", row.names = FALSE)

# --- diseases export ---
# Each disease entry looks like: list(NCIt = "C7961; ...", ORDO = "Orphanet_98834; ...")
split_id_desc <- function(s) {
  s <- trimws(s)
  if (!grepl(";", s, fixed = TRUE)) {
    list(id = s, description = NA_character_)
  } else {
    id  <- sub(";.*$", "", s)
    desc <- sub("^[^;]*;\\s*", "", s)
    list(id = trimws(id), description = trimws(desc))
  }
}

dis_dt <- rbindlist(
  lapply(seq_len(nrow(dt)), function(i) {
    d <- dt$diseases[[i]]
    if (is.null(d) || length(d) == 0) return(NULL)
    acc <- dt$accession[i]
    rbindlist(lapply(names(d), function(src) {
      vals <- d[[src]]
      if (is.null(vals) || length(vals) == 0) return(NULL)
      vals <- as.character(unlist(vals, use.names = FALSE))
      parsed <- lapply(vals, split_id_desc)
      data.table(
        accession = acc,
        source = tolower(src),
        id = vapply(parsed, `[[`, character(1), "id"),
        description = vapply(parsed, `[[`, character(1), "description")
      )
    }), use.names = TRUE, fill = TRUE)
  }),
  use.names = TRUE, fill = TRUE
)

if (!is.null(dis_dt) && nrow(dis_dt) > 0L) {
  # Deduplicate, basic cleanup
  dis_dt[, `:=`(
    accession = trimws(accession),
    source = trimws(source),
    id = trimws(id),
    description = trimws(description)
  )]
  dis_dt <- dis_dt[
    !is.na(accession) & nzchar(accession) &
    !is.na(source) & nzchar(source) &
    !is.na(id) & nzchar(id)
  ]
  dis_dt <- unique(dis_dt, by = c("accession", "source", "id"))
} else {
  dis_dt <- data.table(accession = character(0), source = character(0), id = character(0), description = character(0))
}

# match DB FK/PK column name
if ("accession" %in% names(dis_dt)) {
  data.table::setnames(dis_dt, "accession", "cellosaurus_accession")
}
write.csv(as.data.frame(dis_dt), "output_data/cell_line_diseases.csv", row.names = FALSE)
# --- end diseases export ---

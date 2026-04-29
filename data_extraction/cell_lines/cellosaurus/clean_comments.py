import json
import re
import pandas as pd


def parse_comments(val: str) -> dict:
    """Parse the 'comments' JSON-ish cell safely to a dict."""
    if pd.isna(val) or str(val).strip() in {"", "NA", "NaN", "None"}:
        return {}
    s = str(val)
    try:
        return json.loads(s)
    except Exception:
        try:
            return json.loads(s.replace('""', '"'))
        except Exception:
            return {}


def flatten_value(v):
    """Flatten nested lists/dicts to compact strings suitable for CSV columns."""
    if isinstance(v, list):
        if all(isinstance(x, list) for x in v):
            return " || ".join(" | ".join(map(str, x)) for x in v)
        return " | ".join(map(str, v))
    if isinstance(v, dict):
        return json.dumps(v, ensure_ascii=False)
    if v is None:
        return ""
    return str(v)


SQL_RESERVED_WORDS = {
    "select",
    "from",
    "where",
    "group",
    "order",
    "by",
    "limit",
    "offset",
    "join",
    "on",
    "and",
    "or",
    "not",
    "as",
    "in",
    "is",
    "table",
    "column",
    "index",
    "primary",
    "key",
    "unique",
    "null",
    "true",
    "false",
    "create",
    "drop",
    "insert",
    "update",
    "delete",
}

POST_NORMALIZE_OVERRIDES = {
    "celllinename": "cell_line_name",
    "ageatsampling": "age_at_sampling",
    "sexofcell": "sex_of_cell",
    "crossreferences": "cross_references",
    "hlatyping": "hla_typing",
    "genomeancestry": "genome_ancestry",
    "microsatelliteinstability": "microsatellite_instability",
    "sequencevariation": "sequence_variation",
    "partof": "part_of",
    "karyotypicinformation": "karyotypic_information",
    "problematiccellline": "problematic_cell_line",
    "knockoutcell": "knockout_cell",
    "selectedforresistanceto": "selected_for_resistance_to",
    "group": "group_info",
    "from": "from_col",
}


def _camel_to_snake(name: str) -> str:
    s = re.sub(r"(.)([A-Z][a-z]+)", r"\1_\2", name)
    s = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s)
    return s


def normalize_column_name(col: str) -> str:
    """
    Make a column name SQL-friendly:
      - convert camelCase -> snake_case
      - replace non-alphanum with underscores
      - lowercase, strip leading/trailing underscores
      - prefix underscore if starts with digit
      - avoid reserved words by appending '_col'
      - apply overrides for nicer names
    """
    original = str(col).strip()

    s = _camel_to_snake(original)
    s = re.sub(r"[^\w]+", "_", s)
    s = s.lower().strip("_")

    if not s:
        s = "col"

    if s[0].isdigit():
        s = f"_{s}"

    if s in SQL_RESERVED_WORDS:
        s = f"{s}_col"

    s = POST_NORMALIZE_OVERRIDES.get(s, s)
    return s


def make_unique(names):
    """Ensure column name list is unique by suffixing _2, _3, ... on collisions."""
    seen = {}
    out = []
    for n in names:
        base = n
        i = 2
        while n in seen:
            n = f"{base}_{i}"
            i += 1
        seen[n] = True
        out.append(n)
    return out


def main(in_path: str, out_path: str):
    df = pd.read_csv(in_path)

    idx_cols = [c for c in df.columns if c == "" or str(c).startswith("Unnamed")]
    if idx_cols:
        df = df.drop(columns=idx_cols)

    if "comments" in df.columns:
        comments_series = df["comments"].apply(parse_comments)
        comments_df = comments_series.apply(
            lambda d: {k: flatten_value(v) for k, v in (d or {}).items()}
        ).apply(pd.Series)
        out_df = pd.concat([df.drop(columns=["comments"]), comments_df], axis=1)
    else:
        out_df = df

    normalized = [normalize_column_name(c) for c in out_df.columns]
    normalized = make_unique(normalized)
    out_df.columns = normalized

    out_df.to_csv(out_path, index=False)


if __name__ == "__main__":
    main("output_data/cell_line_table.csv", "output_data/cell_lines_table_cleaned.csv")

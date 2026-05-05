#!/usr/bin/env bash
# Header-aware CSV summarization helpers used by 09_report_generation.sh.
# All functions echo a single token (numeric or string) and never mutate state.
# Sourced by 09_report_generation.sh and tests/unit/test_report_summary.sh.

# csv_col_idx <csv> <col_name>
# Echo the 1-based column index of <col_name>, or "" if not found / file absent.
csv_col_idx() {
    local file="$1" col="$2"
    [ -f "$file" ] || { echo ""; return 0; }
    awk -F',' -v col="$col" 'NR==1{for(i=1;i<=NF;i++)if($i==col){print i;exit}}' "$file"
}

# count_by_value <csv> <col_name> <value>
# Count of data rows where the named column == value (string equality, exact).
count_by_value() {
    local file="$1" col="$2" val="$3"
    local idx
    idx="$(csv_col_idx "$file" "$col")"
    [ -z "$idx" ] && { echo "0"; return 0; }
    awk -F',' -v idx="$idx" -v v="$val" 'NR>1 && $idx==v {n++} END{print n+0}' "$file"
}

# sum_column <csv> <col_name>
# Numeric sum over data rows. Blanks treated as 0; non-numeric tokens contribute 0.
sum_column() {
    local file="$1" col="$2"
    local idx
    idx="$(csv_col_idx "$file" "$col")"
    [ -z "$idx" ] && { echo "0"; return 0; }
    awk -F',' -v idx="$idx" 'NR>1 {s += ($idx + 0)} END{printf "%g", s+0}' "$file"
}

# read_csv_value <csv> <col_name> [row=2]
# Return the value at (row, col). Default row 2 = first data row.
read_csv_value() {
    local file="$1" col="$2" row="${3:-2}"
    local idx
    idx="$(csv_col_idx "$file" "$col")"
    [ -z "$idx" ] && { echo ""; return 0; }
    awk -F',' -v idx="$idx" -v row="$row" 'NR==row {print $idx}' "$file"
}

#!/usr/bin/env bash
# Unit tests for stage 09 report-summary helpers.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
FIXTURE="$PROJECT_ROOT/tests/unit/fixtures/ranked_candidates_may2026_schema.csv"

# shellcheck source=/dev/null
source "$PROJECT_ROOT/scripts/_report_summary_lib.sh"

fail=0
assert_eq() {
    if [ "$1" != "$2" ]; then
        printf 'FAIL: %s — expected %q, got %q\n' "$3" "$2" "$1" >&2
        fail=1
    else
        printf 'PASS: %s\n' "$3"
    fi
}

# csv_col_idx
assert_eq "$(csv_col_idx "$FIXTURE" id)"                  "1"  "csv_col_idx id"
assert_eq "$(csv_col_idx "$FIXTURE" confidence_tier)"     "3"  "csv_col_idx confidence_tier"
assert_eq "$(csv_col_idx "$FIXTURE" selection_significant)" "8" "csv_col_idx selection_significant"
assert_eq "$(csv_col_idx "$FIXTURE" cds_source)"          "35" "csv_col_idx cds_source"
assert_eq "$(csv_col_idx "$FIXTURE" hcr_probe_friendly)"  "39" "csv_col_idx hcr_probe_friendly"
assert_eq "$(csv_col_idx "$FIXTURE" no_such_column)"      ""   "csv_col_idx missing column -> empty"

# count_by_value
assert_eq "$(count_by_value "$FIXTURE" confidence_tier High)"        "1" "count High"
assert_eq "$(count_by_value "$FIXTURE" confidence_tier Medium)"      "2" "count Medium"
assert_eq "$(count_by_value "$FIXTURE" confidence_tier Low)"         "1" "count Low"
assert_eq "$(count_by_value "$FIXTURE" selection_significant True)"  "1" "count selection_significant True"
assert_eq "$(count_by_value "$FIXTURE" busted_s_significant True)"   "1" "count busted_s_significant True"
assert_eq "$(count_by_value "$FIXTURE" busted_mh_significant True)"  "2" "count busted_mh_significant True"
assert_eq "$(count_by_value "$FIXTURE" hcr_probe_friendly True)"     "2" "count hcr_probe_friendly True"
assert_eq "$(count_by_value "$FIXTURE" has_tandem_cluster_data True)" "1" "count has_tandem_cluster_data True"
assert_eq "$(count_by_value "$FIXTURE" cds_source native)"           "2" "count cds_source native"
assert_eq "$(count_by_value "$FIXTURE" cds_source miniprot)"         "1" "count cds_source miniprot"

# sum_column
assert_eq "$(sum_column "$FIXTURE" meme_n_episodic_sites)"           "5" "sum meme_n_episodic_sites"

# read_csv_value
assert_eq "$(read_csv_value "$FIXTURE" id 2)"                "BersteEVm0001t1" "read_csv_value id row 2"
assert_eq "$(read_csv_value "$FIXTURE" confidence_tier 4)"   "Low"             "read_csv_value confidence_tier row 4"

# Missing-file safety: should echo empty / 0, not error.
assert_eq "$(csv_col_idx "/no/such/file" id)"      "" "missing file -> empty idx"
assert_eq "$(count_by_value "/no/such/file" id x)" "0" "missing file -> 0 count"
assert_eq "$(sum_column "/no/such/file" id)"       "0" "missing file -> 0 sum"

if [ "$fail" -eq 0 ]; then
    echo "ALL PASS"
fi
exit $fail

#!/usr/bin/env bash
# Pre-flight validator for the user-supplied external-tool paths referenced
# by an ipg params YAML file. Checks only the licensed / user-provided
# tools — anything bioconda ships is out of scope.
#
# Usage:
#   bin/check_external_tools.sh assets/params_D106_liver.yaml
set -eu

params_file=${1:?"usage: $0 <params.yaml>"}
[[ -r "$params_file" ]] || { echo "cannot read $params_file" >&2; exit 2; }

fail=0
report() {
    local label=$1 path=$2 test=$3
    if [[ -z "$path" || "$path" == "null" ]]; then
        echo "  - $label: unset (skipped)"
    elif [[ "$test" == "executable" ]]; then
        if [[ -x "$path" ]]; then
            echo "  ✓ $label: $path"
        else
            echo "  ✗ $label: $path (missing or not executable)" >&2
            fail=1
        fi
    elif [[ "$test" == "readable" ]]; then
        if [[ -r "$path" ]]; then
            echo "  ✓ $label: $path"
        else
            echo "  ✗ $label: $path (missing or unreadable)" >&2
            fail=1
        fi
    elif [[ "$test" == "blastdb" ]]; then
        if [[ -r "${path}.pin" && -r "${path}.phr" && -r "${path}.psq" ]]; then
            echo "  ✓ $label: ${path}.{pin,phr,psq}"
        else
            echo "  ✗ $label: ${path} (missing .pin/.phr/.psq sibling)" >&2
            fail=1
        fi
    fi
}

yaml_get() {
    # Minimal yaml scalar extractor. Returns empty string if key absent.
    awk -v key="$1" '
        BEGIN { FS=":" }
        $1 == key {
            sub(/^[^:]*:[[:space:]]*/, "")
            gsub(/^[\x27"]|[\x27"]$/, "")
            print; exit
        }
    ' "$params_file"
}

echo "External-tool check for: $params_file"
report "MSFragger JAR"       "$(yaml_get msfragger_jar)"     readable
report "netMHCpan binary"    "$(yaml_get netmhcpan_path)"    executable
report "netMHCIIpan binary"  "$(yaml_get netmhciipan_path)"  executable
report "GibbsCluster script" "$(yaml_get gibbscluster_path)" readable
report "BLAST DB prefix"     "$(yaml_get blast_db)"          blastdb
report "PEAKS PSM CSV"       "$(yaml_get peaks_psm_csv)"     readable

if (( fail )); then
    echo
    echo "One or more required tools are missing. Fix the paths above before kicking off a run." >&2
    exit 1
fi
echo
echo "All configured external tools are present and readable."

#!/usr/bin/env bash
# Fix mzMLs that MSFragger rejects.
#
# Problems addressed:
#   1. cvParam self-closing tags without value="" → Batmass-IO NPE
#   2. Missing <componentList> in instrumentConfiguration → MSFragger
#      can't determine analyzer type (ITMS/FTMS) and rejects file
#   3. indexedmzML byte offsets go stale after edits → strip index
#
# Usage:
#   bash bin/fix_mzml_for_msfragger.sh input.mzML [output.mzML]
set -euo pipefail

input="${1:?usage: $0 input.mzML [output.mzML]}"
output="${2:-$input}"

if [[ "$input" != "$output" ]]; then
    cp "$input" "$output"
fi

python3 -c "
import sys, re

with open('$output', 'r') as f:
    lines = f.readlines()

patched_cv = 0
patched_ic = False
with open('$output', 'w') as f:
    in_index = False
    for i, line in enumerate(lines):
        stripped = line.strip()

        # --- Strip indexedmzML wrapper + stale index ---
        if stripped.startswith('<indexedmzML') or stripped == '</indexedmzML>':
            continue
        if stripped.startswith('<indexList'):
            in_index = True
            continue
        if in_index:
            if '</indexList>' in stripped:
                in_index = False
            continue
        if stripped.startswith('<indexListOffset') or stripped.startswith('<fileChecksum'):
            continue

        # --- Fix cvParam tags missing value= ---
        if '<cvParam' in line and '/>' in line and 'value=' not in line:
            line = line.replace('/>', 'value=\"\" />')
            patched_cv += 1

        # --- Inject componentList if instrumentConfiguration has none ---
        # Detect closing </instrumentConfiguration> and check if the
        # preceding block lacked a <componentList>.  If so, inject a
        # generic Orbitrap FTMS configuration that MSFragger recognises.
        if '</instrumentConfiguration>' in stripped and not patched_ic:
            # Look back up to 10 lines for <componentList>
            lookback = ''.join(lines[max(0,i-10):i])
            if '<componentList' not in lookback:
                indent = line[:len(line) - len(line.lstrip())]
                component = (
                    f'{indent}  <componentList count=\"3\">\\n'
                    f'{indent}    <source order=\"1\">\\n'
                    f'{indent}      <cvParam cvRef=\"MS\" accession=\"MS:1000073\" name=\"electrospray ionization\" value=\"\"/>\\n'
                    f'{indent}    </source>\\n'
                    f'{indent}    <analyzer order=\"2\">\\n'
                    f'{indent}      <cvParam cvRef=\"MS\" accession=\"MS:1000484\" name=\"orbitrap\" value=\"\"/>\\n'
                    f'{indent}    </analyzer>\\n'
                    f'{indent}    <detector order=\"3\">\\n'
                    f'{indent}      <cvParam cvRef=\"MS\" accession=\"MS:1000624\" name=\"inductive detector\" value=\"\"/>\\n'
                    f'{indent}    </detector>\\n'
                    f'{indent}  </componentList>\\n'
                )
                f.write(component)
                patched_ic = True

        f.write(line)

fixes = []
if patched_cv:
    fixes.append(f'{patched_cv} cvParams')
if patched_ic:
    fixes.append('injected componentList')
print(f'[fix_mzml] {\";\".join(fixes) or \"no changes\"} in $output', file=sys.stderr)
"

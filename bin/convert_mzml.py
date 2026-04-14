#!/usr/bin/env python3
"""Convert mzML file(s) to MGF + build scans / index2scan pickles.

Extracted from immunopeptidomics/core.py read_mzML() + write_MGF() (Scull 2021).
Run one mzML per invocation; the MS_SEARCH subworkflow fans out over files and
merges pickles downstream.
"""
from __future__ import annotations

import argparse
import os
import pickle
import re
import sys
from typing import Tuple

import pymzml


def parse_precursor(spectrum) -> Tuple[float, str]:
    try:
        prec = spectrum.selected_precursors[0]
        return prec["mz"], prec["charge"]
    except Exception:
        s = spectrum.to_string().decode("utf-8")
        mz_m = re.search(r'name="selected ion m/z" value="([\d.]+)"', s)
        z_m = re.search(r'name="charge state" value="(\d+)"', s)
        if not (mz_m and z_m):
            raise RuntimeError(f"Unparseable precursor in spectrum {spectrum.ID}")
        return float(mz_m.group(1)), z_m.group(1)


def convert(mzml_path: str, run: str) -> tuple[dict, dict]:
    scans: dict = {}
    index2scan: dict = {}

    with open(f"{run}.mgf", "w") as mgf:
        reader = pymzml.run.Reader(mzml_path)
        for spectrum in reader:
            prec_mz, charge = parse_precursor(spectrum)
            scan_id = str(spectrum.id_dict.get("scan", spectrum.ID))
            rt_sec = spectrum.scan_time_in_minutes() * 60

            mgf.write(
                "BEGIN IONS\n"
                f"TITLE={run}_scan={scan_id}_z={charge}\n"
                f"RTINSECONDS={rt_sec:.6f}\n"
                f"PEPMASS={prec_mz:.4f}\n"
                f"CHARGE={charge}+\n"
            )
            for mz, intensity in zip(spectrum.mz, spectrum.i):
                mgf.write(f"{mz:.6f}\t{intensity}\n")
            mgf.write("END IONS\n")

            if spectrum["ms level"] != 2:
                continue
            scan = str(spectrum.ID)
            specid = f"{run}_{scan}"
            ion_mobility = spectrum.get("MS:1002815", 0)
            index2scan[f"{run}_{spectrum.index + 1}"] = scan
            scans[specid] = {
                "run": run,
                "scan": scan,
                "spectrum_index": spectrum.index,
                "rt": str(spectrum.scan_time_in_minutes()),
                "mz": str(prec_mz),
                "z": str(charge),
                "im": str(ion_mobility),
                "rescore": f"{run}_scan={scan}_z={charge}_seq=",
            }

    return scans, index2scan


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mzml", help="input mzML file")
    p.add_argument("--run", default=None, help="run basename (default: stem of mzML)")
    args = p.parse_args()

    run = args.run or os.path.splitext(os.path.basename(args.mzml))[0]
    scans, index2scan = convert(args.mzml, run)

    with open(f"{run}.scans.pkl", "wb") as f:
        pickle.dump(scans, f)
    with open(f"{run}.index2scan.pkl", "wb") as f:
        pickle.dump(index2scan, f)

    print(f"[convert_mzml] {run}: {len(scans)} MS2 scans -> {run}.mgf", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Generate a synthetic mzML + matching FASTA for the cryptic-peptide positive control.

Emits theoretical b/y-ion MS2 spectra for one planted CRYPTIC peptide plus a deterministic
set of background peptides, and a FASTA containing exactly those sequences. Search engines
(Comet/Sage) should identify every spectrum trivially (perfect fragment matches); the
background gives Mokapot/FDR enough target+decoy PSMs to fit.

Run in the ms2rescore pixi env (needs psims):
    pixi run -e ms2rescore python tests/data/spike/make_synth_spectra.py \
        --out-mzml tests/data/spike/synth.mzML --out-fasta tests/data/spike/synth_db.fasta

Reproducible: fixed RNG seed, so committed fixtures are regenerable byte-for-byte.
"""
from __future__ import annotations

import argparse
import os
import random

from psims.mzml import MzMLWriter

# Monoisotopic residue masses (Da). Cys is bare (no fixed carbamidomethyl).
RESIDUE = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "C": 103.00919, "L": 113.08406, "I": 113.08406, "N": 114.04293,
    "D": 115.02694, "Q": 128.05858, "K": 128.09496, "E": 129.04259, "M": 131.04049,
    "H": 137.05891, "F": 147.06841, "R": 156.10111, "Y": 163.06333, "W": 186.07931,
}
WATER = 18.0105646863
PROTON = 1.0072764666

# The planted cryptic peptide — every letter is a valid amino acid, so it spells a marker
# that is unmistakable in the final identifications and trivially DNA-encodable for P2.
# 11-mer (not a 9-mer): b1/y1 always fall below Sage's fragment_min_mz=200, so a 9-mer
# loses 2 of its 16 peaks and drops under Sage's min_peaks=15 floor. 11 residues give 20
# b/y peaks -> 18 survive the 200 m/z cut, clearing the floor in both Comet and Sage.
CRYPTIC = "CRYPTICALLY"


def peptide_mass(seq: str) -> float:
    return sum(RESIDUE[a] for a in seq) + WATER


def fragment_peaks(seq: str) -> list[tuple[float, float]]:
    """Singly-charged b and y ions, uniform intensity (centroided)."""
    peaks: list[tuple[float, float]] = []
    n = len(seq)
    for i in range(1, n):
        b = sum(RESIDUE[a] for a in seq[:i]) + PROTON
        y = sum(RESIDUE[a] for a in seq[n - i:]) + WATER + PROTON
        peaks.append((b, 1000.0))
        peaks.append((y, 1000.0))
    peaks.sort()
    return peaks


def precursor_mz(seq: str, charge: int) -> float:
    return (peptide_mass(seq) + charge * PROTON) / charge


def random_peptide(rng: random.Random) -> str:
    # Realistic-ish AA pool. Length >= 10: a 10-mer gives 18 b/y peaks, and b1/y1 fall
    # below Sage's fragment_min_mz=200, so 16 survive the cut and clear min_peaks:15.
    # (9-mers leave only 14 after the cut and Sage silently drops them.)
    pool = "ACDEFGHIKLMNPQRSTVWY"
    length = rng.randint(10, 12)
    while True:
        pep = "".join(rng.choice(pool) for _ in range(length))
        if pep != CRYPTIC:
            return pep


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-mzml", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-samplesheet", default=None,
                    help="ms_search samplesheet to emit (default: samplesheet_ms.csv next "
                         "to --out-mzml). nf-schema validates ms_file against launchDir, so "
                         "the path is written absolute to stay portable across run dirs.")
    ap.add_argument("--out-peaks", default=None,
                    help="synthetic PEAKS db.psms.csv to emit (default: db.psms.csv next to "
                         "--out-mzml). Exercises the --ms_engines peaks branch (CONVERT_PEAKS "
                         "ingests a PEAKS export — the commercial tool can't be run here).")
    ap.add_argument("--n-background", type=int, default=300)
    ap.add_argument("--n-noise", type=int, default=150,
                    help="random-peak spectra so some matches land on decoys, "
                         "seeding the Mokapot/FDR decoy distribution")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)

    # Build the peptide list: planted cryptic first, then unique background peptides.
    peptides = [CRYPTIC]
    seen = {CRYPTIC}
    while len(peptides) < args.n_background + 1:
        pep = random_peptide(rng)
        if pep not in seen:
            seen.add(pep)
            peptides.append(pep)

    # Assign charges (immunopeptides are mostly 2+, some 3+).
    charged = [(pep, 2 if i % 4 else 3) for i, pep in enumerate(peptides)]

    # ---- mzML ----------------------------------------------------------------
    with MzMLWriter(open(args.out_mzml, "wb")) as w:
        w.controlled_vocabularies()
        # A standards-complete preamble. MSFragger's Batmass-IO parser is strict:
        # it rejects a <run> without a defaultInstrumentConfigurationRef ("Scans = 0,
        # corrupted file"), whereas Comet/Sage are lenient. So declare an instrument
        # configuration and reference it from the run.
        w.file_description(["MSn spectrum", "centroid spectrum"], None)
        w.software_list([{"id": "synth_gen", "version": "1.0",
                          "params": ["custom unreleased software tool"]}])
        w.instrument_configuration_list([
            w.InstrumentConfiguration(
                id="IC1",
                component_list=[w.Source(1, ["electrospray ionization"]),
                                w.Analyzer(2, ["orbitrap"]),
                                w.Detector(3, ["inductive detector"])])
        ])
        w.data_processing_list([
            w.DataProcessing(
                processing_methods=[w.ProcessingMethod(
                    order=1, software_reference="synth_gen",
                    params=["Conversion to mzML"])],
                id="DP1")
        ])
        with w.run(id="synth_spike", instrument_configuration="IC1"):
            with w.spectrum_list(count=len(charged) + args.n_noise):
                # Real spectra — clean b/y ions; high-scoring target matches.
                for i, (pep, z) in enumerate(charged, start=1):
                    peaks = fragment_peaks(pep)
                    w.write_spectrum(
                        [p[0] for p in peaks], [p[1] for p in peaks],
                        id=f"scan={i}",
                        centroided=True,
                        scan_start_time=float(i),
                        params=[{"ms level": 2}, "MSn spectrum"],
                        precursor_information={"mz": precursor_mz(pep, z), "charge": z},
                    )
                # Decoy-seed spectra — a FEW real b/y ions of the REVERSED (decoy)
                # sequence (so they match the rev_ decoy at >= min_matched_peaks:4 and
                # within the 10 ppm fragment tolerance) PLUS random filler peaks to
                # clear Sage's min_peaks:15 floor (filler never matches at 10 ppm).
                # The partial match scores BELOW the full-b/y targets, so Mokapot gets
                # a low-scoring decoy distribution and the planted peptide stays on top.
                bg = peptides[1:]  # exclude the planted cryptic peptide
                for j in range(args.n_noise):
                    rev = rng.choice(bg)[::-1]
                    z = 2 if j % 4 else 3
                    real = fragment_peaks(rev)
                    matched = rng.sample(real, 6)                       # 6 true decoy ions
                    filler = [(rng.uniform(100.0, 1200.0), 300.0)       # >=15 total peaks
                              for _ in range(12)]
                    peaks = sorted(matched + filler)
                    scan = len(charged) + j + 1
                    w.write_spectrum(
                        [p[0] for p in peaks], [p[1] for p in peaks],
                        id=f"scan={scan}",
                        centroided=True,
                        scan_start_time=float(scan),
                        params=[{"ms level": 2}, "MSn spectrum"],
                        precursor_information={"mz": precursor_mz(rev, z), "charge": z},
                    )

    # ---- FASTA ---------------------------------------------------------------
    # Each peptide is its own entry. The cryptic marker is tagged so the assertion
    # and origins step can find it; background entries mimic canonical proteins.
    with open(args.out_fasta, "w") as f:
        for i, (pep, _z) in enumerate(charged):
            if pep == CRYPTIC:
                f.write(f">CRYPTIC_SPIKE_{pep}\n{pep}\n")
            else:
                f.write(f">sp|SYN{i:04d}|SYN{i:04d}_HUMAN synthetic background\n{pep}\n")

    # ---- samplesheet ---------------------------------------------------------
    # nf-schema resolves the ms_file 'file-path' field against launchDir, not the
    # samplesheet dir, so write an absolute path (regenerated per machine).
    samplesheet = args.out_samplesheet or os.path.join(
        os.path.dirname(os.path.abspath(args.out_mzml)), "samplesheet_ms.csv")
    with open(samplesheet, "w") as f:
        f.write("sample,ms_file,condition,fraction,replicate\n")
        f.write(f"SPIKE,{os.path.abspath(args.out_mzml)},synthetic,1,rep1\n")

    # ---- PEAKS db.psms.csv (exercises the --ms_engines peaks branch) ----------
    # PEAKS is commercial and can't be run here; the pipeline ingests its db.psms.csv
    # export (CONVERT_PEAKS -> mokapot). Emit a synthetic export: target PSMs (scans 1..N,
    # incl. CRYPTICALLY at scan 1 with a high -10LgP) + decoy PSMs (reversed background
    # peptides, low -10LgP) so Mokapot/Percolator can fit FDR. Scan = mzML scan index
    # (index2scan from CONVERT_MZML maps synth_<scan> -> <scan>).
    peaks_csv = args.out_peaks or os.path.join(
        os.path.dirname(os.path.abspath(args.out_mzml)), "db.psms.csv")
    src = os.path.basename(args.out_mzml)
    bg = peptides[1:]
    with open(peaks_csv, "w") as f:
        f.write("Source File,Scan,decoy,Accession,-10LgP,Tag Length,Delta RT,"
                "MS2 Correlation,ppm,Peptide\n")
        for i, (pep, _z) in enumerate(charged, start=1):              # target PSMs
            acc = f"CRYPTIC_SPIKE_{pep}" if pep == CRYPTIC else f"sp|SYN{i - 1:04d}|SYN{i - 1:04d}_HUMAN"
            lgp = round(rng.uniform(45, 70) if pep == CRYPTIC else rng.uniform(30, 60), 2)
            f.write(f"{src},{i},False,{acc},{lgp},{len(pep)},{round(rng.uniform(-0.5, 0.5), 3)},"
                    f"{round(rng.uniform(0.85, 0.98), 3)},{round(rng.uniform(-2, 2), 2)},{pep}\n")
        for j in range(args.n_noise):                                 # decoy PSMs (low scores)
            scan = len(charged) + j + 1
            rev = rng.choice(bg)[::-1]
            f.write(f"{src},{scan},True,#DECOY#rev_SYN{j:04d},{round(rng.uniform(8, 28), 2)},"
                    f"{len(rev)},{round(rng.uniform(-3, 3), 3)},{round(rng.uniform(0.2, 0.45), 3)},"
                    f"{round(rng.uniform(-6, 6), 2)},{rev}\n")

    print(f"[make_synth_spectra] wrote {len(charged) + args.n_noise} spectra "
          f"({len(charged)} target + {args.n_noise} decoy-seed) to {args.out_mzml}")
    print(f"[make_synth_spectra] wrote {len(charged)} sequences to {args.out_fasta}")
    print(f"[make_synth_spectra] wrote samplesheet to {samplesheet}")
    print(f"[make_synth_spectra] wrote PEAKS db.psms.csv to {peaks_csv} "
          f"({len(charged)} target + {args.n_noise} decoy PSMs)")
    print(f"[make_synth_spectra] planted cryptic peptide: {CRYPTIC}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

Changes made in Revision 1 of the research note. The revised note text is with the journal and is not in this repository; the figure, new runs, scripts and source files it cites are in this folder. Paths in this file are the author's working paths.

# AAS82162 rev1 (TIC 71150324): changes from submission5/rnaas5.md

Decision: REVISE. With no neighbours removed, FPP = 0.018 ± 0.003 (10/10 runs ≤ 0.035), which meets the 0.035 bar with no exclusions at all.

## New TRICERATOPS + MOLUSC runs (revision-82162/fpp/)
- Script `fpp/run.py`, copied from `prerev/71150324/fpp/run.py` (the submitted joint run). Frozen: TRICERATOPS 1.1.0 (`triceratops/.venv`), pickle `deep/71150324/fpp/out/target_71150324.pkl` (45 stars), Gaia background, MOLUSC file `molusc_for_trc.csv`, filt TESS, N = 10⁶, exptime 30/1440, nsamples 20, seeds 1000–1009, 10-min phase bins, vet duration 0.0817 d.
- Changed: the light curve is all six sectors (`deep/71150324/fit/all6/fitdata.npz`: S5, 32, 98, 106, 107, 108; 1,631 points, wotan biweight with a transit mask). It is folded on the fitted P = 8.422169 d and T0 = 4256.4492 BTJD, and tdepth is the fitted 2277.4 ppm (was 1787.4 from S5/32/98 only). Phase window ±0.15 d, the fit-data extent (44 bins); err = 838 ppm.
- `fpp/summary_full.json`. No neighbours removed: FPP 0.0185 ± 0.0033 (0.0125–0.0232), NFPP 0.0062 ± 0.0011. Scenarios: TP 0.839, PTP 0.141, STP 0.0117, NEB 0.0058, DTP 0.0013.
- `fpp/perstar_nfpp_full.csv`. NFPP shares: TIC 71150320 94.4%, 671498495 3.7%, 671498530 1.9%, 71150328 0.002%.
- Pixel rule (summed Δχ² ≥ 25 and no sector preferring the neighbour by ≥ 25), from `prerev/71150324/loc/dchi2_sum_6sectors{,_persector}.csv`. TIC 71150320: Σ 112.2, min sector +3.0. 671498530: Σ 95.2, min +2.7. 671498495: Σ 77.6, min +2.0. All three meet it.
- `fpp/summary_sel.json`. Only those 3 removed (`drop_sel.txt`): FPP 0.0123 ± 0.0031 (0.0070–0.0168), NFPP 2.7 × 10⁻⁷ (all from TIC 71150328, 84.8″, no pixel test).
- Note: full_0–7 ran 4 at a time (`go.sh`); full_8–9 and sel_0–9 ran 6 at a time (`go2.sh`) after the first runner hit its time limit. Same script and seeds.

## Text changes
1. Address line removed (name, "Independent researcher", email only).
2. Abstract: the FPP is now the no-removal value, 0.018 ± 0.003, "above the 0.015 validation threshold". Localisation is "all six sectors". "Same-CCD comparison stars" and "four independent extractions" were cut for length.
3. Search: the lightkurve Savitzky–Golay flatten (0.75 d, no mask, 30-min bins; `hunt7/run1.py` lines 36–38) is now stated, and the search is separated from the fit. The search sectors are S5, 32, 98. The three search gates were cut.
4. Fit: wotan biweight 0.5 d is named, with ±2.9 h masked and the trend interpolated (`fit/all6/prep.py` lines 22–23; mask = 3 × 0.0817 d total width).
5. Vetting: the per-transit depth and error are now defined (`sectors/scripts/measure.py` events/comb). The √χ²ᵣ inflation is per pooled set (`systematics/inflated_errors.json`: 8.42σ, o/e 0.065σ, phase 0.5 −69 ± 388), and the text now says so. The comparison-star fractions are given as numbers (`systematics/sametime_summary.json`, `s106_108_comparisons.json`). The momentum-dump sentence is now "within 0.1 d" (the computed quantity; accepted S106–108 transits are ≥ 0.58 d from dumps, `new_sector_spoc_dumps.json`). The single-transit median S/N of 3.6 is from `prerev/71150324/fig/fig_stats.json`. The raw-vs-corrected and post-egress sentences were cut for length.
6. Localisation: three sectors became six. Δχ² is defined (free PRF amplitude + constant at each position, noise from 120 random-phase difference images; `hunt7/localize/71150324/diff.py` lines 37–51). Six-sector sums are 77.6–112.4, single sectors 1.2–35.3, and S108 is weak (≤ 3.9).
7. FPP section rewritten: the no-removal result comes first, then the NFPP carriers, the pixel justification, and the removal of only those three (Giacalone et al. 2021). MOLUSC use is spelled out (200,000 generated, 50.8% kept; `prerev/71150324/methods_molusc_counts.json`). "Gaia astrometry and imaging" became "Gaia RUWE and contrast tests". The text states that sub-0.015 is reached only with pixel exclusions and that there is no imaging, so no validation claim.
8. Figure 1: new three panels (folded; unfolded S98 QLP with trend; four single transits, S5/32/107/108), from `prerev/71150324/fig/`. The caption is rewritten.
9. Reference added: Hippke et al. 2019 (wotan).

Abstract 141 words; body (Star → Summary) 899 words.

## Not changed / open
- S109 is in the footprint but not in the data (no light curve in `sectors/`).

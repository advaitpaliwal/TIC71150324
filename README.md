# TIC 71150324.01 — a temperate 2.5 Earth-radius planet candidate on an 8.4-day orbit

Discovered by **Advait Paliwal** (independent researcher, advaitpaliwal@gmail.com) on 2026-10-10 in TESS full-frame-image photometry.

| | |
|---|---|
| Host | TIC 71150324 = 2MASS J04293004−1654557, M dwarf, 95.0 pc, T = 12.80 |
| Period | 8.4221690 (+0.0000171/−0.0000158) d |
| Mid-transit | 2461256.4492 ± 0.0027 BJD_TDB |
| Depth / duration | 2277 ± 300 ppm / 1.84 h |
| Radius | 2.49 ± 0.19 R⊕ |
| Equilibrium temperature | 480 ± 21 K |
| TESS data | sectors 5, 32, 98, 106, 107, 108 (FFI only); 17 transits |
| False-positive probability | 0.018 ± 0.003 (TRICERATOPS + MOLUSC Gaia limits, no neighbours removed, fitted depth; Revision 1); NFPP 0.0062 |

Status: **planet candidate**. Not a TOI or CTOI as of 2026-10-10. High-resolution imaging is needed to exclude an unresolved bound companion.

## Revision 1 (2026-10-16)

The research note was revised in response to the journal editor. The corrected numbers are:

- **False-positive probability, no neighbours removed:** FPP = 0.0185 ± 0.0033 (ten runs, range 0.0125–0.0232), NFPP = 0.0062 ± 0.0011, all 45 stars kept (`revision1/fpp/summary_full.json`). TRICERATOPS was given the fitted depth of 2277 ppm (was 1787 ppm) and the six-sector fit light curve folded on the fitted ephemeris. This is above the 0.015 validation threshold. The version 1.0 value of 0.033 ± 0.003 used pixel exclusions and is replaced.
- **NFPP carriers:** TIC 71150320 (94%), TIC 671498495 (4%) and TIC 671498530 (2%). Each is disfavoured by the six-sector difference images with summed Δχ² ≥ 77.6 and no sector preferring it.
- **Removing only those three** (following Giacalone et al. 2021): FPP = 0.0123 ± 0.0031 (range 0.0070–0.0168), NFPP = 2.7 × 10⁻⁷ (`revision1/fpp/summary_sel.json`). This rests on pixel exclusions and there is no imaging, so no validation is claimed.
- **Difference images:** six sectors (was three); summed Δχ² 77.6–112.4.
- **MOLUSC:** Gaia RUWE and contrast tests only (no imaging); 200,000 companions generated, 50.8% kept.
- **Figure 1** now has three panels: folded, unfolded Sector 98 with the trend, and four single transits (`revision1/figure1_combined_TIC71150324.png`).

Every change, with its source, is listed in [`revision1/CHANGES.md`](revision1/CHANGES.md). The revised note text is with the journal and is not included here.

Data: TESS (NASA, via MAST: TESS-SPOC doi:10.17909/t9-wpz1-8s54, QLP doi:10.17909/t9-r086-e880) and ESA Gaia. Licences: code MIT; text, figures and data CC BY 4.0.

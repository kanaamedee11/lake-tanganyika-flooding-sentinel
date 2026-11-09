# Where the Lake Rises: Sentinel-1 and Sentinel-2 evidence of shoreline and delta flooding at northern Lake Tanganyika, 2017 to 2025

Lake Tanganyika reached record levels in 2020 to 2021 and again in 2024. What does free Sentinel imagery show at the lake's northern end
(Bujumbura, Gatumba, the Rusizi delta)? *Status: preprint draft, not peer reviewed.* Paper: [`paper/paper.pdf`](paper/paper.pdf).

## Main findings
* **Lake-wide monthly water area is nearly uninformative.** Sentinel-1 and Sentinel-2 differ by 22 km2 RMSD on a lake of about 340 km2 and are uncorrelated month to month, more than the 9 km2 rise in persistent water between 2019 and 2025.
* **The change is local.** A low-lying area about 2 km south-west of Gatumba went from water in 1% of radar passes in 2017 to 92% in 2025. Sentinel-2 independently shows the same ramp in the same pixels (2% in 2019 to 99% in 2025, correlation 0.995 with the radar), and it persists at other radar thresholds and in both orbit directions.
* The connected area is 4.5 km2 (radar) and 6.2 km2 (Sentinel-2). Newly persistent water totals about 7 km2 (4 to 10 km2 across threshold and period choices).
* **Limits:** no field or lake-level data, one site, and sensors that detect open water but not flooded buildings.

## Layout
`lake_tanganyika_export.js`, `lake_tanganyika_export2.js` (Earth Engine exports) -> `data/` -> `analysis.py`, `analysis2.py` -> `results/` -> `make_figs.py`, `make_figs2.py` -> `paper/`.

## Reproduce
1. `pip install -r requirements.txt`.
2. (Optional) re-run the two Earth Engine scripts in the Code Editor and put the exports in `data/`. The exports used in the paper are already included.
3. From the repository root: `python analysis.py`, `python analysis2.py`, `python make_figs.py`, `python make_figs2.py`.
4. `cd paper && pdflatex paper && bibtex paper && pdflatex paper && pdflatex paper`.

Thresholds are fixed (Sentinel-1 VV below -16 dB, MNDWI above 0); see the limitations section of the paper. Data terms: `DATA.md`.

## Licence and citation
Code MIT (`LICENSE`); paper text and figures CC BY 4.0 (`LICENSE-PAPER.md`); cite via `CITATION.cff`. Contact: Kana Amédée, kanaamedee00@gmail.com.

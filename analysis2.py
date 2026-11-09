"""Second round of checks: Sentinel-2 pixel maps, threshold sensitivity, orbit split, west extension."""
import json, numpy as np, rasterio
from scipy import ndimage as ndi, stats
D = "data/"; OUT = "results/"
PX = 900 / 1e6; YEARS = list(range(2017, 2026))
def rd(name):
    a = rasterio.open(D + name).read().astype(float); a[a == 255] = np.nan; return a
s1 = rd("tanganyika_s1_water_frequency_by_year.tif"); s2 = rd("tanganyika_s2_frequency_by_year.tif")
thr = rd("tanganyika_s1_threshold_sensitivity.tif"); orb = rd("tanganyika_s1_by_orbit_direction.tif")
occ = rasterio.open(D + "tanganyika_reference_layers.tif").read()[0]
M = np.load(OUT + "masks.npz"); big = M["big"]; gain16 = M["gain"]
R = {}
def mean_in(a, m): return [float(np.nanmean(a[i][m])) if np.isfinite(a[i][m]).any() else None for i in range(a.shape[0])]
def gainmask(a, bi, li, cut=50):
    b = np.nanmean(a[bi], axis=0); l = np.nanmean(a[li], axis=0); v = np.isfinite(b) & np.isfinite(l)
    return v & (l >= cut) & (b < cut), b, l
def biggest(g):
    lab, n = ndi.label(g); 
    if n == 0: return g & False, 0.0
    sz = ndi.sum(g, lab, range(1, n + 1)); k = np.argmax(sz) + 1; return lab == k, float(sz.max() * PX)

# ---- 1. Sentinel-2 confirmation in the main cluster
R["s2_nodata_share_by_year"] = [float(np.isnan(s2[i]).mean()) for i in range(9)]
R["s2_cluster_freq"] = mean_in(s2, big); R["s1_cluster_freq"] = mean_in(s1, big)
R["s2_cluster_valid_share"] = [float(np.isfinite(s2[i][big]).mean()) for i in range(9)]
ok = [i for i in range(9) if R["s2_cluster_freq"][i] is not None]
a = np.array([R["s1_cluster_freq"][i] for i in ok]); b = np.array([R["s2_cluster_freq"][i] for i in ok])
R["cluster_s1_s2_pearson_years"] = float(stats.pearsonr(a, b)[0]); R["cluster_years_used"] = [YEARS[i] for i in ok]
# S2-only change map: base 2019-2020, late 2024-2025 (2017-18 have too few clear observations)
g2, b2, l2 = gainmask(s2, slice(2, 4), slice(7, 9))
R["s2_gain_km2"] = float(g2.sum() * PX); big2, R["s2_main_km2"] = biggest(g2)
# fair S1 comparison over the same periods
g1b, _, _ = gainmask(s1, slice(2, 4), slice(7, 9)); R["s1_gain_same_periods_km2"] = float(g1b.sum() * PX)
valid = np.isfinite(b2) & np.isfinite(l2) & np.isfinite(np.nanmean(s1[2:4], axis=0)) & np.isfinite(np.nanmean(s1[7:9], axis=0))
inter = (g1b & g2 & valid).sum(); union = ((g1b | g2) & valid).sum()
R["s1_s2_gain_iou"] = float(inter / union); R["s1_s2_gain_overlap_of_s2_share"] = float(inter / max(1, (g2 & valid).sum()))
R["s2_gain_inside_s1_cluster_share"] = float((g2 & big).sum() / max(1, g2.sum()))
# S2 frequency in the cluster by year restricted to 2019 on, plus the S1 curve at same years
# ---- 2. S1 threshold sensitivity (cluster mask from the -16 dB analysis stays fixed; gain recomputed per threshold)
t14 = thr[:9]; t18 = thr[9:]
R["thr_cluster_freq"] = {"-14": mean_in(t14, big), "-16": mean_in(s1, big), "-18": mean_in(t18, big)}
for name, a_ in (("-14", t14), ("-16", s1), ("-18", t18)):
    g, _, _ = gainmask(a_, slice(0, 3), slice(6, 9)); bg, km = biggest(g)
    R.setdefault("thr_gain", {})[name] = {"gain_km2": float(g.sum() * PX), "main_km2": km, "persist_2019": float(np.nansum(a_[2] >= 50) * PX), "persist_2025": float(np.nansum(a_[8] >= 50) * PX)}
deep = (occ >= 90) & (ndi.distance_transform_edt(occ >= 90) * 30 > 1000)
R["thr_deep_freq"] = {"-14": mean_in(t14, deep), "-16": mean_in(s1, deep), "-18": mean_in(t18, deep)}
# ---- 3. by orbit direction
asc = orb[:9]; desc = orb[9:]
R["orbit_cluster_freq"] = {"asc": mean_in(asc, big), "desc": mean_in(desc, big)}
R["orbit_deep_freq"] = {"asc": mean_in(asc, deep), "desc": mean_in(desc, deep)}
R["orbit_nodata_share"] = {"asc": [float(np.isnan(asc[i]).mean()) for i in range(9)], "desc": [float(np.isnan(desc[i]).mean()) for i in range(9)]}
for name, a_ in (("asc", asc), ("desc", desc)):
    g, _, _ = gainmask(a_, slice(0, 3), slice(6, 9)); bg, km = biggest(g); R.setdefault("orbit_gain", {})[name] = {"gain_km2": float(g.sum() * PX), "main_km2": km}
# ---- 4. west extension (grid-aligned with the main box: x offset 740 px west, y offset 367 px north)
w = rasterio.open(D + "tanganyika_delta_west_extension.tif").read().astype(float); w[w == 255] = np.nan
w1, w2 = w[:9], w[9:]
for name, a_ in (("s1", w1), ("s2", w2)):
    if name == "s1": g, _, _ = gainmask(a_, slice(0, 3), slice(6, 9))
    else: g, _, _ = gainmask(a_, slice(2, 4), slice(7, 9))
    bg, km = biggest(g)
    R.setdefault("west", {})[name] = {"gain_km2": float(g.sum() * PX), "main_km2": km, "main_freq_by_year": mean_in(a_, bg)}
    if name == "s1": west_big = bg
# the connected area in the extended box vs the old truncated area: how far west does it go and how big is it overall
ys, xs = np.where(west_big)
R["west"]["s1_main_cols_px"] = [int(xs.min()), int(xs.max())]; R["west"]["s1_main_touches_west_edge"] = bool(xs.min() == 0)
R["west"]["s1_main_touches_east_edge"] = bool(xs.max() == w1.shape[2] - 1)
# S1 persistent area in extended box by year
R["west"]["s1_persistent_km2"] = [float(np.nansum(w1[i] >= 50) * PX) for i in range(9)]
json.dump(R, open(OUT + "results2.json", "w"), indent=1, default=float)
print(json.dumps(R, indent=1, default=float))

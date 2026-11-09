"""Lake Tanganyika (Bujumbura / Rusizi delta) water-extent analysis from the user's Earth Engine exports."""
import sys, json, numpy as np, pandas as pd, rasterio
from scipy import stats, ndimage as ndi
from rasterio.warp import transform
DATA = "data/"
OUT = "results/"
PX = 900 / 1e6  # km2 per 30 m pixel
R = {}

df = pd.read_csv(DATA + "tanganyika_monthly_water_timeseries.csv")
df["date"] = pd.to_datetime(df.month + "-01"); df["year"] = df.date.dt.year
FULL = 1500.0   # km2 of clear (valid) pixels required, about 88% of the 1,700 km2 AOI
df["s2ok"] = df.s2_valid_km2 >= FULL; df["jrcok"] = df.jrc_valid_km2 >= FULL
R["n_months"] = len(df); R["n_s2ok"] = int(df.s2ok.sum()); R["n_jrcok"] = int(df.jrcok.sum())
R["aoi_km2"] = float(df.jrc_valid_km2.max())

def agree(a, b, mask):
    d = df.loc[mask, [a, b]].dropna(); x, y = d[a].values, d[b].values
    return {"n": int(len(d)), "bias": float((x - y).mean()), "rmsd": float(np.sqrt(((x - y) ** 2).mean())),
            "mad": float(np.abs(x - y).mean()), "pearson": float(stats.pearsonr(x, y)[0]),
            "spearman": float(stats.spearmanr(x, y)[0]), "mean_b": float(y.mean())}
R["s1_vs_s2"] = agree("s1_water_km2", "s2_water_km2", df.s2ok)
R["s1_vs_jrc"] = agree("s1_water_km2", "jrc_water_km2", df.jrcok)
R["s2_vs_jrc"] = agree("s2_water_km2", "jrc_water_km2", df.s2ok & df.jrcok)
d = df[df.s2ok].copy(); d["diff"] = d.s1_water_km2 - d.s2_water_km2
R["s1_minus_s2_quantiles"] = {k: float(v) for k, v in d["diff"].quantile([.05, .25, .5, .75, .95]).items()}
R["worst_s1_low_months"] = d.nsmallest(5, "diff")[["month", "s1_water_km2", "s2_water_km2", "diff", "s1_scenes"]].round(1).to_dict("records")
R["s1_all_low_months"] = df.nsmallest(5, "s1_water_km2")[["month", "s1_water_km2", "s1_scenes"]].round(1).to_dict("records")
# do dips in S1 coincide with S2 flat? months where S1 is >=20 km2 below its own median
med = df.s1_water_km2.median(); R["s1_median"] = float(med)
R["s1_dips_gt20"] = int((df.s1_water_km2 < med - 20).sum())

# trends: Theil-Sen on monthly series with seasonal-month anomalies, moving-block bootstrap CI
def theil_block(t, y, nboot=2000, block=12, seed=0):
    ok = ~np.isnan(y); t, y = t[ok], y[ok]
    sl = stats.theilslopes(y, t)[0]
    res = y - (stats.theilslopes(y, t)[1] + sl * t)
    rng = np.random.default_rng(seed); n = len(y); bs = []
    for _ in range(nboot):
        starts = rng.integers(0, max(1, n - block), size=int(np.ceil(n / block)))
        idx = np.concatenate([np.arange(s, min(s + block, n)) for s in starts])[:n]
        yy = (stats.theilslopes(y, t)[1] + sl * t) + res[idx]
        bs.append(stats.theilslopes(yy, t)[0])
    return float(sl), [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))], int(len(y))
df["ty"] = (df.date.dt.year - 2017) + (df.date.dt.month - 1) / 12
def deseason(col, mask):
    s = df[col].where(mask); clim = s.groupby(df.date.dt.month).transform("mean"); return s - clim
R["trend_s1_km2_per_yr"] = theil_block(df.ty.values, deseason("s1_water_km2", df.s1_scenes > 0).values)
R["trend_s2_km2_per_yr"] = theil_block(df.ty.values, deseason("s2_water_km2", df.s2ok).values)
yr = df.groupby("year").agg(s1=("s1_water_km2", "median"), s1_scenes=("s1_scenes", "sum"), n_s2ok=("s2ok", "sum"))
yr["s2"] = df[df.s2ok].groupby("year").s2_water_km2.median()
R["annual_median"] = yr.round(1).reset_index().to_dict("records")

# rasters
ds = rasterio.open(DATA + "tanganyika_s1_water_frequency_by_year.tif"); tr = ds.transform
occ, lc, el, sl_ = rasterio.open(DATA + "tanganyika_reference_layers.tif").read()
fr = ds.read().astype(float); excl = (fr == 255); fr[excl] = np.nan
R["excluded_share"] = float(excl[0].mean())
years = list(range(2017, 2026))
persist = {y: float(np.nansum(fr[i] >= 50) * PX) for i, y in enumerate(years)}
inter = {y: float(np.nansum((fr[i] >= 10) & (fr[i] < 50)) * PX) for i, y in enumerate(years)}
R["persistent_km2"] = persist; R["intermittent_km2"] = inter
R["jrc_perm_km2"] = float((occ >= 50).sum() * PX)

def gainmask(b_idx, l_idx, cut=50):
    base = np.nanmean(fr[b_idx], axis=0); late = np.nanmean(fr[l_idx], axis=0)
    v = ~np.isnan(base) & ~np.isnan(late)
    return v & (late >= cut) & (base < cut), v & (base >= cut) & (late < cut), base, late
g, l, base, late = gainmask(slice(0, 3), slice(6, 9))
R["gain_km2"] = float(g.sum() * PX); R["loss_km2"] = float(l.sum() * PX)
lab, n = ndi.label(g); sizes = ndi.sum(g, lab, range(1, n + 1)); big = lab == (np.argmax(sizes) + 1)
R["n_gain_clusters"] = int(n); R["main_cluster_km2"] = float(big.sum() * PX)
ys, xs = np.where(big); cx, cy = tr * (xs.mean(), ys.mean()); lo, la = transform(ds.crs, "EPSG:4326", [cx], [cy])
R["main_cluster_lonlat"] = [float(lo[0]), float(la[0])]
R["main_cluster_elev"] = np.percentile(el[big], [5, 50, 95]).tolist()
R["main_cluster_freq_by_year"] = {y: float(np.nanmean(fr[i][big])) for i, y in enumerate(years)}
rest = g & ~big
R["rest_gain_km2"] = float(rest.sum() * PX); R["rest_freq_by_year"] = {y: float(np.nanmean(fr[i][rest])) for i, y in enumerate(years)}
R["gain_landcover_km2"] = {int(k): float(v * PX) for k, v in zip(*np.unique(lc[g], return_counts=True))}
R["main_landcover_km2"] = {int(k): float(v * PX) for k, v in zip(*np.unique(lc[big], return_counts=True))}
dist = ndi.distance_transform_edt(~(occ >= 90)) * 30
R["gain_dist_perm_m"] = np.percentile(dist[g], [10, 50, 90]).tolist()
# sensitivity of the gain estimate to cut and period choice
sens = []
for name, bi, li in (("2017-19 vs 2023-25", slice(0, 3), slice(6, 9)), ("2017-18 vs 2024-25", slice(0, 2), slice(7, 9)),
                     ("2017-19 vs 2023-24", slice(0, 3), slice(6, 8)), ("2017-19 vs 2022-23", slice(0, 3), slice(5, 7))):
    for cut in (30, 50, 70):
        gg, ll, _, _ = gainmask(bi, li, cut)
        lb, nn = ndi.label(gg); ss = ndi.sum(gg, lb, range(1, nn + 1)) if nn else [0]
        sens.append({"periods": name, "cut": cut, "gain": float(gg.sum() * PX), "loss": float(ll.sum() * PX), "main": float(max(ss) * PX)})
R["sensitivity"] = sens
# agreement of S1 persistent water with JRC long-term occurrence (not time matched)
mid = np.nanmean(fr[2:5], axis=0); v = ~np.isnan(mid); pred = mid >= 50; truth = occ >= 50
tp = (pred & truth & v).sum(); fp = (pred & ~truth & v).sum(); fn = (~pred & truth & v).sum()
R["s1_vs_jrc_occ"] = {"iou": float(tp / (tp + fp + fn)), "precision": float(tp / (tp + fp)), "recall": float(tp / (tp + fn)),
                      "fp_km2": float(fp * PX), "fn_km2": float(fn * PX)}
# same but for the earliest years (before most of the rise)
mid0 = np.nanmean(fr[0:3], axis=0); v0 = ~np.isnan(mid0); p0 = mid0 >= 50
tp = (p0 & truth & v0).sum(); fp = (p0 & ~truth & v0).sum(); fn = (~p0 & truth & v0).sum()
R["s1_vs_jrc_occ_2017_19"] = {"iou": float(tp / (tp + fp + fn))}
json.dump(R, open(OUT + "results.json", "w"), indent=1, default=float)
np.savez_compressed(OUT + "masks.npz", gain=g, loss=l, big=big, base=base, late=late)
print(json.dumps({k: R[k] for k in ["n_months", "n_s2ok", "n_jrcok", "s1_vs_s2", "s1_vs_jrc", "s2_vs_jrc", "s1_minus_s2_quantiles", "worst_s1_low_months", "s1_dips_gt20", "trend_s1_km2_per_yr", "trend_s2_km2_per_yr"]}, indent=1, default=float))
print(json.dumps({k: R[k] for k in ["annual_median", "persistent_km2", "gain_km2", "loss_km2", "main_cluster_km2", "main_cluster_lonlat", "gain_dist_perm_m", "s1_vs_jrc_occ", "s1_vs_jrc_occ_2017_19", "excluded_share"]}, indent=1, default=float))
for s in sens: print(s)

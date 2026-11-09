import json, numpy as np, pandas as pd, rasterio
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
HALO=[pe.withStroke(linewidth=2, foreground='white')]
from matplotlib.colors import LightSource, ListedColormap, LinearSegmentedColormap
from matplotlib.patches import Rectangle
from rasterio.warp import transform
DATA = "data/"; OUT = "paper/figures/"
R = json.load(open("results/results.json")); M = np.load("results/masks.npz")
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"], "font.size": 8,
                     "axes.labelsize": 8, "axes.titlesize": 8, "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7, "pdf.fonttype": 42})
ds = rasterio.open(DATA + "tanganyika_s1_water_frequency_by_year.tif"); tr = ds.transform
occ, lc, el, sl = rasterio.open(DATA + "tanganyika_reference_layers.tif").read()
fr = ds.read().astype(float); fr[fr == 255] = np.nan
H, W = occ.shape
def lonlat(col, row):
    x, y = tr * (col, row); lo, la = transform(ds.crs, "EPSG:4326", [x], [y]); return lo[0], la[0]
l0, t0 = lonlat(0, 0); l1, b1 = lonlat(W, H)
EXT = [l0, l1, b1, t0]
def px_of(lon, lat):
    x, y = transform("EPSG:4326", ds.crs, [lon], [lat]); c, r = ~tr * (x[0], y[0]); return c, r
BUJ = (29.36, -3.38); GAT = (29.25, -3.333)
ls = LightSource(azdeg=315, altdeg=40); elc = np.clip(el.astype(float), 765, 1700)
hs = ls.hillshade(elc, vert_exag=3, dx=30, dy=30)
ys, xs = np.where(M["big"]); pad = 45
c0, c1, r0, r1 = max(xs.min() - pad, 0), min(xs.max() + pad, W), max(ys.min() - pad, 0), min(ys.max() + pad, H)

# ---------- Fig 1 study area
fig, ax = plt.subplots(figsize=(3.45, 3.7))
ax.imshow(hs, cmap="gray", extent=EXT, vmin=0.2, vmax=1.0)
ov = np.ma.masked_where(occ < 50, np.ones_like(occ)); ax.imshow(ov, cmap=ListedColormap(["#4a90d9"]), extent=EXT, alpha=0.85)
for (lo, la), name, dx, dy in ((BUJ, "Bujumbura", 0.01, 0.0), (GAT, "Gatumba", 0.012, 0.01)):
    ax.plot(lo, la, "o", ms=3.5, mfc="yellow", mec="k", mew=0.5); ax.text(lo + dx, la + dy, name, fontsize=7, color="black", weight="bold", path_effects=HALO)
a, b = lonlat(c0, r0); c, d = lonlat(c1, r1)
ax.add_patch(Rectangle((a, d), c - a, b - d, fill=False, ec="red", lw=1.0)); ax.text(a, d - 0.018, "Fig. 3 and 4 area", color="red", fontsize=6.5, path_effects=HALO)
ax.text(29.27, -3.50, "Lake Tanganyika", color="white", fontsize=8, style="italic")
ax.set_xlabel("Longitude (E)"); ax.set_ylabel("Latitude"); ax.set_xlim(EXT[0], EXT[1]); ax.set_ylim(EXT[2], EXT[3])
plt.tight_layout(); plt.savefig(OUT + "fig_area.pdf", bbox_inches="tight"); plt.savefig(OUT + "fig_area.png", dpi=250, bbox_inches="tight"); plt.close()

# ---------- Fig 2 monthly series + scatter
df = pd.read_csv(DATA + "tanganyika_monthly_water_timeseries.csv"); df["date"] = pd.to_datetime(df.month + "-01")
s2ok = df.s2_valid_km2 >= 1500; jok = df.jrc_valid_km2 >= 1500
fig = plt.figure(figsize=(7.0, 2.7)); g = fig.add_gridspec(1, 2, width_ratios=[2.2, 1], wspace=0.28)
a1 = fig.add_subplot(g[0]); a2 = fig.add_subplot(g[1])
a1.plot(df.date, df.s1_water_km2, "-", color="#555555", lw=0.8, label="Sentinel-1 (all months)")
a1.plot(df.date[s2ok], df.s2_water_km2[s2ok], "o", ms=2.6, color="#0072B2", label="Sentinel-2 (clear months)")
a1.plot(df.date[jok], df.jrc_water_km2[jok], "s", ms=2.6, mfc="none", mec="#D55E00", mew=0.7, label="JRC monthly (clear months)")
a1.axvline(pd.Timestamp("2021-12-23"), color="k", ls=":", lw=0.7); a1.text(pd.Timestamp("2021-12-23"), 366, " S1B stops", fontsize=6.5, va="top")
a1.set_ylim(275, 380); a1.set_ylabel("Water area in study box (km$^2$)"); a1.set_title("(a) Monthly lake-wide water area"); a1.grid(alpha=0.3, lw=0.4)
a1.legend(loc="upper left", ncol=3, frameon=True, framealpha=0.9, fontsize=6.3, columnspacing=0.8, handletextpad=0.3)
m = s2ok
a2.plot(df.s2_water_km2[m], df.s1_water_km2[m], "o", ms=2.6, color="#0072B2", alpha=0.8); a2.plot([275, 365], [275, 365], "k--", lw=0.6)
a2.set_xlim(325, 365); a2.set_ylim(275, 365); a2.set_xlabel("Sentinel-2 (km$^2$)"); a2.set_ylabel("Sentinel-1 (km$^2$)"); a2.set_title("(b) Same month, two sensors"); a2.grid(alpha=0.3, lw=0.4)
plt.savefig(OUT + "fig_series.pdf", bbox_inches="tight"); plt.savefig(OUT + "fig_series.png", dpi=250, bbox_inches="tight"); plt.close()

# ---------- Fig 3 small multiples around the main cluster
yrs = [2017, 2019, 2021, 2023, 2025]; cm = LinearSegmentedColormap.from_list("w", ["#f4f1e8", "#9ecae1", "#08306b"])
fig, axs = plt.subplots(1, 5, figsize=(7.0, 1.5), constrained_layout=True)
lo0, la0 = lonlat(c0, r0); lo1, la1 = lonlat(c1, r1)
for ax, y in zip(axs, yrs):
    w = fr[y - 2017][r0:r1, c0:c1]; ax.set_facecolor("#c9c9c9")
    im = ax.imshow(np.ma.masked_invalid(w), cmap=cm, vmin=0, vmax=100, extent=[lo0, lo1, la1, la0]); ax.set_title(str(y), fontsize=8)
    ax.set_xticks([]); ax.set_yticks([])
    ax.plot(*GAT, "o", ms=3, mfc="yellow", mec="k", mew=0.4) if lo0 < GAT[0] < lo1 and la1 < GAT[1] < la0 else None
cb = fig.colorbar(im, ax=axs, fraction=0.02, pad=0.01); cb.set_label("Passes classed as water (%)", fontsize=7); cb.ax.tick_params(labelsize=6.5)
plt.savefig(OUT + "fig_delta_years.pdf", bbox_inches="tight"); plt.savefig(OUT + "fig_delta_years.png", dpi=250, bbox_inches="tight"); plt.close()

# ---------- Fig 4 gain/loss map + curves
fig = plt.figure(figsize=(7.0, 2.9)); g = fig.add_gridspec(1, 2, width_ratios=[1.15, 1], wspace=0.25)
a1 = fig.add_subplot(g[0]); a2 = fig.add_subplot(g[1])
sub = (slice(r0, r1), slice(c0, c1)); EX = [lo0, lo1, la1, la0]
a1.imshow(hs[sub], cmap="gray", extent=EX, vmin=0.3, vmax=1.0)
a1.imshow(np.ma.masked_where(occ[sub] < 50, np.ones_like(occ[sub])), cmap=ListedColormap(["#9ecae1"]), extent=EX, alpha=0.9)
a1.imshow(np.ma.masked_where(~M["gain"][sub], np.ones(M["gain"][sub].shape)), cmap=ListedColormap(["#d62828"]), extent=EX)
a1.imshow(np.ma.masked_where(~M["loss"][sub], np.ones(M["loss"][sub].shape)), cmap=ListedColormap(["#f4a261"]), extent=EX)
a1.plot(*GAT, "o", ms=4, mfc="yellow", mec="k", mew=0.5); a1.text(GAT[0] - 0.012, GAT[1] - 0.004, "Gatumba", fontsize=7, color="black", weight="bold", path_effects=HALO)
a1.set_xlabel("Longitude (E)"); a1.set_ylabel("Latitude"); a1.set_title("(a) New persistent water, 2017-19 to 2023-25")
from matplotlib.lines import Line2D
a1.legend(handles=[Line2D([], [], marker="s", ls="", color="#d62828", label="gained water"), Line2D([], [], marker="s", ls="", color="#f4a261", label="lost water"),
                   Line2D([], [], marker="s", ls="", color="#9ecae1", label="JRC permanent water")], loc="lower left", fontsize=6.5, framealpha=0.9)
yrs_all = list(range(2017, 2026))
a2.plot(yrs_all, [R["main_cluster_freq_by_year"][str(y)] for y in yrs_all], "-o", ms=3, color="#d62828", label="main cluster (4.5 km$^2$)")
a2.plot(yrs_all, [R["rest_freq_by_year"][str(y)] for y in yrs_all], "-s", ms=3, color="#555555", label="other gained pixels (2.5 km$^2$)")
for y, t in ((2020.4, "2020/21 high water"), ): pass
a2.axvspan(2019.5, 2021.5, color="#0072B2", alpha=0.12); a2.axvspan(2023.5, 2024.5, color="#0072B2", alpha=0.12)
a2.text(2020.5, 4, "reported\nfloods\n2020-21", ha="center", fontsize=6.5, color="#0072B2"); a2.text(2024.0, 4, "reported\nrecord\n2024", ha="center", fontsize=6.5, color="#0072B2")
a2.set_xlabel("Year"); a2.set_ylabel("Passes classed as water (%)"); a2.set_ylim(0, 100); a2.set_title("(b) Mean yearly water frequency"); a2.grid(alpha=0.3, lw=0.4)
a2.legend(loc="upper left", fontsize=6.5)
plt.savefig(OUT + "fig_gain.pdf", bbox_inches="tight"); plt.savefig(OUT + "fig_gain.png", dpi=250, bbox_inches="tight"); plt.close()
print("ok", EXT, (c0, c1, r0, r1))

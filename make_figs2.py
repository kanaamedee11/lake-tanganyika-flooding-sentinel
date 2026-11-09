import json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R = json.load(open("results/results2.json")); OUT = "paper/figures/"
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"], "font.size": 8,
                     "axes.labelsize": 8, "axes.titlesize": 8, "legend.fontsize": 6.8, "xtick.labelsize": 7, "ytick.labelsize": 7, "pdf.fonttype": 42})
Y = list(range(2017, 2026))
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.7), sharey=True)
for k, c, ls, lab in (("-14", "#56B4E9", "--", "radar, $-14$ dB"), ("-16", "#0072B2", "-", "radar, $-16$ dB (main)"), ("-18", "#999999", ":", "radar, $-18$ dB")):
    a1.plot(Y, R["thr_cluster_freq"][k], ls, color=c, marker="o", ms=2.5, label=lab)
s2 = [np.nan if v is None else v for v in R["s2_cluster_freq"]]
a1.plot(Y, s2, "-", color="#D55E00", marker="s", ms=3, lw=1.3, label="Sentinel-2, MNDWI $>0$")
a1.set_title("(a) Two sensors and three radar thresholds"); a1.set_xlabel("Year"); a1.set_ylabel("Observations classed as water (%)")
a1.legend(loc="upper left", frameon=True, framealpha=0.9); a1.grid(alpha=0.3, lw=0.4); a1.set_ylim(0, 102)
a2.plot(Y, R["orbit_cluster_freq"]["asc"], "-", color="#CC79A7", marker="^", ms=3, label="gained area, ascending")
a2.plot(Y, R["orbit_cluster_freq"]["desc"], "-", color="#009E73", marker="v", ms=3, label="gained area, descending")
a2.plot(Y, R["orbit_deep_freq"]["asc"], "--", color="#CC79A7", label="open lake, ascending")
a2.plot(Y, R["orbit_deep_freq"]["desc"], "--", color="#009E73", label="open lake, descending")
a2.set_title("(b) Radar split by orbit direction"); a2.set_xlabel("Year"); a2.legend(loc="center left", frameon=True, framealpha=0.9, bbox_to_anchor=(0.0, 0.62)); a2.grid(alpha=0.3, lw=0.4)
plt.tight_layout(w_pad=1.0); plt.savefig(OUT + "fig_checks.pdf", bbox_inches="tight"); plt.savefig(OUT + "fig_checks.png", dpi=250, bbox_inches="tight")

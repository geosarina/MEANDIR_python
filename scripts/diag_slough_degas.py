"""Investigate the slough degassing discrepancy: per-sim distribution of the
degassing contribution to DIC for the 4 slough samples, scenario 2.
"""
import numpy as np
from meandir.run import run

UE = "reference/data/MEANDIR_UserEntries.xlsx"
RD = "reference/data/RiverDataSpreadsheet_Kemeny_etal_2023.xlsx"
SLOUGH = [14, 16, 17, 23]

r, ctx = run("AK_scenario2_carbonate_slctindi_degas_2p5", UE, RD, seed=1,
             max_success=60, max_iter=8000, max_zerohits=4000, samples=SLOUGH)
p = ctx["params"]
di = p.EMList0.index("degas")

print("Per-slough-sample degassing contribution to DIC (%):")
for s in SLOUGH:
    succ = r.sample_results[s]
    if not succ:
        print(f"  sample {s}: 0 successes")
        continue
    dic = np.array([x["fractions"]["DIC"][di] * 100 for x in succ])
    print(f"  sample {s}: n={len(succ)}  median={np.median(dic):+.1f}%  "
          f"[5th {np.percentile(dic, 5):+.1f}, 95th {np.percentile(dic, 95):+.1f}]  "
          f"frac near 0 (>=-0.2%): {100 * np.mean(dic >= -0.2):.0f}%")

alld = np.array([x["fractions"]["DIC"][di] * 100 for s in SLOUGH for x in r.sample_results[s]])
print(f"POOLED slough: median={np.median(alld):+.1f}%  5th={np.percentile(alld, 5):+.1f}  "
      f"95th={np.percentile(alld, 95):+.1f}  band={np.percentile(alld, 95) - np.percentile(alld, 5):.0f}pts "
      f"(pub -3.9)  frac near 0: {100 * np.mean(alld >= -0.2):.0f}%")

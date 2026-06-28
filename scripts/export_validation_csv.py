"""Export validation-test results to CSV under validation_outputs/.

Currently produces the headline Table S2 comparison (fractional contributions,
200 successes/sample) from the per-sample cache written by validate_groups.py.
Run validate_groups first (so the cache exists), then this.
"""
import csv
import json
import os

import numpy as np

from scripts.validate_groups import PUB, classify, aggregate, SCEN

CACHE = "/tmp/claude-0/-home-user-MEANDIR-python/feec57dc-e16b-557c-afd4-47f9242d36db/scratchpad/validate_groups_cache_s200.json"
OUTDIR = "validation_outputs"
DEGAS = {1: "none", 2: "<2.5x DIC", 3: "<25x DIC"}
GROUPKEY = {"All samples": "All", "Mainstem": "Main", "Slough": "Slough"}


def export_tableS2():
    cache = json.load(open(CACHE))
    main, slough = classify()
    groups = {"All samples": main + slough, "Mainstem": main, "Slough": slough}
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, "tableS2_comparison_200successes.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["group", "n_samples", "scenario", "degassing_constraint",
                    "ion", "end_member", "published_pct", "python_pct", "delta_pct"])
        for gname, gidx in groups.items():
            agg = {sc: aggregate(cache, sc, gidx) for sc in SCEN}
            for sc in (1, 2, 3):
                for ion, lab, em, pub in PUB:
                    t = pub[GROUPKEY[gname]][sc - 1]
                    v = agg[sc][(ion, em)]
                    vs = "" if (v != v) else round(v, 2)
                    ts = "" if t is None else t
                    ds = "" if (t is None or v != v) else round(v - t, 2)
                    w.writerow([gname, len(gidx), f"S{sc}", DEGAS[sc],
                                ion, lab, ts, vs, ds])
    print(f"wrote {path} ({sum(1 for _ in open(path)) - 1} data rows)")


if __name__ == "__main__":
    export_tableS2()

#!/usr/bin/env python3
import numpy as np
from scipy import stats
from modules.rank0_scalar_pole import evaluate_rank0
from modules.rank1_foliation_null import evaluate_rank1
from modules.rank2_curvature_dihiggs import evaluate_rank2
from modules.rank3_cartan_triality import evaluate_rank3
from modules.rank4_vbs_unitarization import evaluate_rank4

def main():
    print("=" * 86)
    print("MODULAR RADION (602.50 GeV) MULTITENSOR ANALYSIS SUITE")
    print("=" * 86)

    r0 = evaluate_rank0()
    r1 = evaluate_rank1()
    r2 = evaluate_rank2()
    r3 = evaluate_rank3()
    r4 = evaluate_rank4()

    results = [r0, r1, r2, r3, r4]
    header = "{:<6} | {:<38} | {:<10} | {:<10} | p-value".format("Rank", "Physics Sector", "Delta Chi2", "Z (sigma)")
    print(header)
    print("-" * 86)
    for r in results:
        line = "Rank {:<1} | {:<38} | {:<10.2f} | {:<10.2f} | {:.3e}".format(
            r["rank"], r["name"], r["delta_chi2"], r["z_score"], r["p_val"]
        )
        print(line)
    print("-" * 86)

    delta_chi2_global = r0["delta_chi2"] + r2["delta_chi2"] + r3["delta_chi2"] + r4["delta_chi2"]
    z_global = float(np.sqrt(delta_chi2_global))
    p_global = float(stats.chi2.sf(delta_chi2_global, df=4))

    print("")
    print("GLOBAL COMBINED SIGNAL SIGNIFICANCE (Ranks 0 + 2 + 3 + 4):")
    print("  Total Profile Likelihood Delta Chi2:  {:.2f}".format(delta_chi2_global))
    print("  Combined Statistical Significance:   Z_global = {:.2f} sigma".format(z_global))
    print("  Global Asymptotic p-value:           p_global = {:.2e}".format(p_global))
    print("=" * 86)

if __name__ == "__main__":
    main()

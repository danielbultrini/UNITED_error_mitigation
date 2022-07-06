# %%
import pandas as pd


def relative_shots(training_sets, noise_levels, max_copies, shots_budget):
    CDR_shots = shots_budget // training_sets
    vnCDR_shots = CDR_shots // noise_levels
    CGVD_shots = shots_budget // training_sets
    storage = pd.DataFrame(
        {
            "VD": shots_budget // 2,
            "ZNE": shots_budget // 2,
            "CDR": [CDR_shots] * 6,
            "vnCDR": [vnCDR_shots] * 6,
            "CGVD": [
                base_cost // (copy + 1) // 2
                for copy, base_cost in enumerate([CGVD_shots] * max_copies)
            ],
            "vnCGVD": [
                x // noise_levels // 2
                for x in [
                    base_cost // (copy + 1)
                    for copy, base_cost in enumerate([CGVD_shots] * max_copies)
                ]
            ],
            "copies": list(range(1, max_copies + 1)),
        }
    )
    return pd.melt(storage, id_vars=("copies"))


def get_budgets(shot_budgets=[10**5, 10**6, 10**7, 10**8, 10**9, 10**10]):
    budgets = []
    for i in shot_budgets:
        budgets.extend(relative_shots(100, 3, 6, i)["value"].unique().tolist())
    return sorted(list(set(budgets)))


# %%

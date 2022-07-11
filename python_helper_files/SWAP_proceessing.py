# %%
import pandas as pd
import matplotlib.pyplot as plt
from post_processing_multi import *
import seaborn as sns


def relative_shots(training_sets, noise_levels, max_copies, shots_budget):
    if shots_budget >= 10**9:
        CDR_shots = shots_budget // training_sets
        vnCDR_shots = CDR_shots // noise_levels
        CGVD_shots = shots_budget // training_sets
        storage = pd.DataFrame(
            {
                "VD": [CDR_shots] * 6,
                "ZNE": [CDR_shots] * 6,
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
    elif shots_budget != 0:
        CDR_shots = shots_budget // training_sets
        vnCDR_shots = CDR_shots // noise_levels
        CGVD_shots = shots_budget // training_sets
        storage = pd.DataFrame(
            {
                "VD": [shots_budget // 2] * (max_copies),
                "ZNE": [shots_budget // 2] * (max_copies),
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

    else:
        storage = pd.DataFrame(
            {
                "VD": None,
                "ZNE": None,
                "CDR": [None] * 6,
                "vnCDR": [None] * 6,
                "CGVD": [None] * 6,
                "vnCGVD": [None] * 6,
                "copies": list(range(1, max_copies + 1)),
            }
        )

    return pd.melt(storage, id_vars=("copies"))


def relative_shots_table(training_sets, noise_levels, max_copies, shots_budget):
    if shots_budget != 0:
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
    else:
        storage = pd.DataFrame(
            {
                "VD": None,
                "ZNE": None,
                "CDR": None,
                "vnCDR": None,
                "CGVD": None,
                "vnCGVD": None,
                "copies": list(range(1, max_copies + 1)),
            }
        )

    return storage


def single_qubit_gate_length(layers, qubits):
    return layers * 6 * (qubits - 1)


def two_qubit_gate_length(layers, qubits):
    e2 = layers * (qubits - 1)


def circuit_length(layers, qubits):
    e1 = layers * 6 * (qubits - 1)
    e2 = layers * (qubits - 1)
    return e1 + e2


def total_error_rates_normal(layers, qubits, factor=1):
    e1 = layers * 6 * (qubits - 1) * 0.0011 * factor
    e2 = layers * (qubits - 1) * 0.0021 * factor
    return e1 + e2


def load_multiple_files(
    Q,
    p,
    num_seeds,
    nlsp_list=[1, 2, 3],
    N=10,
    Nts=100,
    max_copies=6,
    tags=[""],
    density_matrices=False,
    shots=[None],
    folder="",
    train=True,
):
    """Special class tha considers swap noise levels"""
    dfs_error_rescaling = []
    dfs_absolute_error = []
    dfs_standard_values = []
    for qubit_no, depth in zip(Q, p):
        print(f"Qubit: {Q}, depth: {p}")
        data_processing = VD_CDR(
            qubit_no,
            depth,
            num_seeds,
            N,
            Nts,
            max_copies,
            nlsp_list,
            shots=shots,
            folder=folder,
            autoprocessing=True,
            plots=False,
            density_matrices=density_matrices,
            extra_tags=tags,
            train=train,
        )

        dfs_absolute_error.append(
            data_processing.abs_error_df.assign(Qubits=qubit_no, depth=depth)
        )
        dfs_error_rescaling.append(
            data_processing.calculated_vals_df.assign(Qubits=qubit_no, depth=depth)
        )
        dfs_standard_values.append(
            data_processing.standard_df.assign(Qubits=qubit_no, depth=depth)
        )
    dfs_absolute_error = pd.concat(dfs_absolute_error)
    dfs_error_rescaling = pd.concat(dfs_error_rescaling)
    dfs_standard_values = pd.concat(dfs_standard_values)

    return dfs_absolute_error, dfs_error_rescaling, dfs_standard_values


def load_multiple_files_budget(
    Q,
    p,
    num_seeds,
    nlsp_list=[1, 2, 3],
    N=10,
    Nts=100,
    max_copies=6,
    tags=[""],
    density_matrices=False,
    shots=[[None]],
    folders=[""],
    train=True,
    budgets=[1],
    train_use=100,
):
    """Special class that considers swap noise levels"""
    dfs_error_rescaling = []
    dfs_absolute_error = []
    dfs_standard_values = []
    for tag, folder in zip(tags, folders):
        print(tag)
        for budget_i, budget in enumerate(budgets):
            for qubit_no, depth in zip(Q, p):
                print(f"Qubit: {qubit_no}, depth: {depth}, budget:{budget}")
                data_processing = VD_CDR(
                    qubit_no,
                    depth,
                    num_seeds,
                    N,
                    Nts,
                    max_copies,
                    nlsp_list,
                    shots=shots[budget_i],
                    folder=folder,
                    autoprocessing=True,
                    density_matrices=density_matrices,
                    extra_tags=[tag],
                    train=train,
                    budget=budget,
                    train_use=train_use,
                )
                dfs_absolute_error.append(
                    data_processing.abs_error_df.assign(
                        Qubits=qubit_no, depth=depth, budget=budget
                    )
                )
                dfs_error_rescaling.append(
                    data_processing.calculated_vals_df.assign(
                        Qubits=qubit_no, depth=depth, budget=budget
                    )
                )
                dfs_standard_values.append(
                    data_processing.standard_df.assign(
                        Qubits=qubit_no, depth=depth, budget=budget
                    )
                )
    dfs_absolute_error = pd.concat(dfs_absolute_error)
    dfs_error_rescaling = pd.concat(dfs_error_rescaling)
    dfs_standard_values = pd.concat(dfs_standard_values)

    return dfs_absolute_error, dfs_error_rescaling, dfs_standard_values


def filter_budget(df, budgets):
    """Special class to filter results with swap noise levels."""
    filtered_dfs = []
    for budget_i, budget in enumerate(budgets):
        shot_budget = relative_shots(100, 3, 6, budget)
        for deets in shot_budget.iterrows():
            copy_no = deets[1][0]
            technique = deets[1][1]
            allowance = deets[1][2]
            if budget == 0:
                df_temp_abs = df.query(
                    f'budget=={budget} & type.str.startswith("{technique}") & copies == {copy_no}'
                )
            else:
                df_temp_abs = df.query(
                    f'shots == {allowance} & type.str.startswith("{technique}") & copies == {copy_no}'
                ).assign(budget=budget)
            filtered_dfs.append(df_temp_abs)

    filtered_dfs = pd.concat(filtered_dfs)
    return filtered_dfs


# %%
base_folder = "./results/"
folders = [
    base_folder
]  # ,base_folder+'LinGrow_complete/realistic_noise_model_shots/']#,base_folder+'RNM_e0p1/',base_folder+'RNM_e1/']
tags = ["_nlswp1"]
budgies = [
    [None],
    # [50000,  1000,   333,   500,   250,   166,   125,   100,    83, 55,    41,    33,    27],
    #  [500000,  10000,   3333,   5000,   2500,   1666,   1250,   1000, 833,    555,    416,    333,    277],
    #  [5000000,  100000,   33333,   50000,   25000,   16666,   12500, 10000,    8333,    5555,    4166,    3333,    2777] ,
    #  [50000000,  1000000,   333333,   500000,   250000,   166666, 125000,   100000,    83333,    55555,    41666,    33333, 27777]
]
budge_abs, budge_res, budge_val = load_multiple_files_budget(
    [4],
    [5],
    30,
    nlsp_list=[1, 2, 3],
    N=30,
    Nts=200,
    max_copies=6,
    tags=tags,
    density_matrices=False,
    shots=budgies,
    folders=folders,
    train=True,
    budgets=[0, 1, 2, 3, 4],
    train_use=200,
)

# %%
t = filter_budget(budge_abs, [0, 10**5, 10**6, 10**7, 10**8, 10**9, 10**10])
t["volume"] = t.depth * t.Qubits * t.Qubits


t.to_pickle(base_folder + "swap_nosie_filtered.pkl")

# t = pd.read_pickle('./BIG_CONGLOMERATE/'+"EVERYTHING_no_swap_nosie_filtered.pkl")
# %%
scaled_SWAP_Results = []
for i in [0, 10**5, 10**6, 10**7, 10**8, 10**9, 10**10]:
    for QUBIT, p in zip([4], [5]):
        temp = budge_abs.query(
            f"abs_error != 0 & nlsp == 1 & budget == {i} & Qubits == {QUBIT} & depth == {p}"
        ).copy()
        temp["scaled_error"] = (
            temp["abs_error"]
            / budge_abs.query(
                f'abs_error != 0 & nlsp == 1 & budget == {i} & copies == 1 & type == "VD_LinGrow" & Qubits == {QUBIT} & depth == {p}'
            )["abs_error"].mean()
        )
        scaled_SWAP_Results.append(temp)
scaled_swap_panda = pd.concat(scaled_SWAP_Results)
scaled_swap_panda.to_pickle("./BIG_CONGLOMERATE/" + "swap_nosie_scaled_filtered.pkl")

# %%
tags = ["_nlswp1"]
title_dict = {
    "_nlswp1": "Noisy swaps",
    "_RNM_e1": "Scaled Noise (e_r=1)",
    "_RNM_e0p1": "Scaled Noise (e_r=0.1)",
}
for tag in tags:
    for QUBIT, g in zip([4], [5]):
        fig, axes = plt.subplots(2, 3, figsize=(15, 10), sharey="row", sharex=True)
        axes = axes.flatten()
        dfs = []
        for i in [10**7, 10**8, 10**9, 10**10]:
            df1 = t.query(
                f"budget=={i}    & Qubits == {QUBIT} & depth == {g} & abs_error != 0 "
            )
            df2 = t.query(
                f"budget=={i*10}    & Qubits == {QUBIT} & depth == {g} & abs_error != 0"
            )
            dfs.append(pd.concat([df1, df2]))
        sns.lineplot(
            data=t.query(
                f'budget==0 &   Qubits == {QUBIT} & depth == {g} & abs_error != 0      &  type.str.endswith("{tag}")'
            ),
            ax=axes[0],
            x="copies",
            y="abs_error",
            hue="type",
            estimator="mean",
            err_style=None,
        )
        sns.lineplot(
            data=t.query(
                f'budget==0 &   Qubits == {QUBIT} & depth == {g} & abs_error != 0      &  type.str.endswith("{tag}")'
            ),
            ax=axes[3],
            x="copies",
            y="abs_error",
            hue="type",
            estimator="mean",
            err_style=None,
        )

        sns.lineplot(
            data=dfs[0].query(
                f'   Qubits == {QUBIT} & depth == {g} & abs_error != 0      &  type.str.endswith("{tag}") '
            ),
            ax=axes[0 + 1],
            x="copies",
            y="abs_error",
            hue="type",
            estimator="mean",
            err_style=None,
        )
        sns.lineplot(
            data=dfs[1].query(
                f'    Qubits == {QUBIT} & depth == {g} & abs_error != 0 &  type.str.endswith("{tag}") '
            ),
            ax=axes[1 + 1],
            x="copies",
            y="abs_error",
            hue="type",
            estimator="mean",
            err_style=None,
        )
        sns.lineplot(
            data=dfs[2].query(
                f' Qubits == {QUBIT} & depth == {g} & abs_error != 0 &  type.str.endswith("{tag}")'
            ),
            ax=axes[2 + 2],
            x="copies",
            y="abs_error",
            hue="type",
            estimator="mean",
            err_style=None,
        )
        sns.lineplot(
            data=dfs[3].query(
                f' Qubits == {QUBIT} & depth == {g} & abs_error != 0 &  type.str.endswith("{tag}") '
            ),
            ax=axes[3 + 2],
            x="copies",
            y="abs_error",
            hue="type",
            estimator="mean",
            err_style=None,
        )
        for i in axes:
            i.set_yscale("log")
            # i.set_ylim([10**(-2),10**1])
            i.legend([], [], frameon=False)

        axes[1].set_title("10^7   Budget")
        axes[2].set_title("10^8   Budget")
        axes[4].set_title("10^9   Budget")
        axes[5].set_title("10^10   Budget")

        handles, labels = axes[0].get_legend_handles_labels()

        fig.legend(
            handles,
            labels,
            title="Mitigation Strategies:",
            loc="lower center",
            bbox_to_anchor=(0.5, -0.025),
            ncol=7,
            fancybox=True,
            shadow=True,
        )

        # axes[2].set_title('1M  Budget')
        # axes[3].set_title('10M Budget')
        fig.suptitle(f"Noisy swap at {title_dict[tag]} Q={QUBIT} depth={g}")

# %%

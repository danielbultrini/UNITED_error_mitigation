# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# plt.style.use(['science','no-latex'])
from post_processing_multi import *
import seaborn as sns

sns.set_context("paper")


def relative_shots(training_sets, noise_levels, max_copies, shots_budget):
    if shots_budget >= 10**9:
        CDR_shots = shots_budget // training_sets
        vnCDR_shots = CDR_shots // noise_levels
        CGVD_shots = shots_budget // training_sets
        storage = pd.DataFrame(
            {
                "VD": [shots_budget] + [shots_budget // 2] * (max_copies - 1),
                "ZNE": [shots_budget // 2] + [shots_budget // 4] * (max_copies - 1),
                "CDR": [CDR_shots] + [CDR_shots // 2] * (max_copies - 1),
                "vnCDR": [vnCDR_shots] * (max_copies),
                "CGVD": [
                    base_cost // (copy + 1) // 2
                    for copy, base_cost in enumerate([CGVD_shots] * (max_copies))
                ],
                "vnCGVD": [
                    x // noise_levels // 2
                    for x in [
                        base_cost // (copy + 1)
                        for copy, base_cost in enumerate([CGVD_shots] * (max_copies))
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
                "ZNE": [shots_budget // 2] + [shots_budget // 4] * (max_copies - 1),
                "CDR": [CDR_shots] + [CDR_shots // 2] * (max_copies - 1),
                "vnCDR": [vnCDR_shots] * (max_copies),
                "CGVD": [
                    base_cost // (copy + 1) // 2
                    for copy, base_cost in enumerate([CGVD_shots] * (max_copies))
                ],
                "vnCGVD": [
                    x // noise_levels // 2
                    for x in [
                        base_cost // (copy + 1)
                        for copy, base_cost in enumerate([CGVD_shots] * (max_copies))
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
                "VD": [shots_budget] + [shots_budget // 2] * (max_copies - 1),
                "ZNE": [shots_budget // 2] + [shots_budget // 4] * (max_copies - 1),
                "CDR": [CDR_shots] + [CDR_shots // 2] * (max_copies - 1),
                "vnCDR": [vnCDR_shots] + [vnCDR_shots // 2] * (max_copies - 1),
                "CGVD": [CGVD_shots]
                + [
                    base_cost // (copy + 2) // 2
                    for copy, base_cost in enumerate([CGVD_shots] * (max_copies - 1))
                ],
                "vnCGVD": [vnCDR_shots]
                + [
                    x // noise_levels // 2
                    for x in [
                        base_cost // (copy + 2)
                        for copy, base_cost in enumerate(
                            [CGVD_shots] * (max_copies - 1)
                        )
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


def cost_fucntion_eval_table(training_sets, noise_levels, max_copies, shots_budget):
    CDR_shots = shots_budget // training_sets
    vnCDR_shots = CDR_shots // noise_levels
    CGVD_shots = shots_budget // training_sets
    storage = pd.DataFrame(
        {
            "VD": [1] + [2] * (max_copies - 1),
            "ZNE": [2] + [4] * (max_copies - 1),
            "CDR": [100] + [200] * (max_copies - 1),
            "vnCDR": [300] + [600] * (max_copies - 1),
            "CGVD": [100]
            + [
                base_cost * (copy + 2) * 2
                for copy, base_cost in enumerate([100] * (max_copies - 1))
            ],
            "vnCGVD": [300]
            + [
                x * noise_levels * 2
                for x in [
                    base_cost * (copy + 2)
                    for copy, base_cost in enumerate([100] * (max_copies - 1))
                ]
            ],
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
    dfs_error_rescaling = []
    dfs_absolute_error = []
    dfs_standard_values = []
    for tag, folder in zip(tags, folders):
        print(tag)
        for budget_i, budget in enumerate(budgets):
            for qubit_no, depth in zip(Q, p):
                print(f"Qubit: {qubit_no}, depth: {depth}, budget:{budget}")
                # try:
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
                # except:
                #     print("FAILED")
                #     pass
    dfs_absolute_error = pd.concat(dfs_absolute_error)
    dfs_error_rescaling = pd.concat(dfs_error_rescaling)
    dfs_standard_values = pd.concat(dfs_standard_values)
    # dfs_absolute_error.to_pickle( folder+f"abs_error.pkl")
    # dfs_error_rescaling.to_pickle(folder+f"error_rescaling.pkl")
    # dfs_standard_values.to_pickle(folder+f"raw_results.pkl")

    return dfs_absolute_error, dfs_error_rescaling, dfs_standard_values


def filter_budget(df, budgets):
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
# budge_abs = pd.read_pickle('./new_data/BIG_CONGLOMERATE_SHOTS/BIG_CONGLOMERATE/'+"abs_error.pkl")

base_folder = "D:/Databases/VD_CDR/SWAP_experiments/"
folders = [
    base_folder + "/swap_data_p_35/"
]  # ,base_folder+'LinGrow_complete/realistic_noise_model_shots/']#,base_folder+'RNM_e0p1/',base_folder+'RNM_e1/']
tags = ["_nlsw1", "_nlsw0"]
budgies = [
    [None],
    [50000, 1000, 333, 500, 250, 166, 125, 100, 83, 55, 41, 33, 27],
    [500000, 10000, 3333, 5000, 2500, 1666, 1250, 1000, 833, 555, 416, 333, 277],
    [
        5000000,
        100000,
        33333,
        50000,
        25000,
        16666,
        12500,
        10000,
        8333,
        5555,
        4166,
        3333,
        2777,
    ],
    [
        50000000,
        1000000,
        333333,
        500000,
        250000,
        166666,
        125000,
        100000,
        83333,
        55555,
        41666,
        33333,
        27777,
    ],
    [
        500000000,
        10000000,
        3333333,
        5000000,
        2500000,
        1666666,
        1250000,
        1000000,
        833333,
        555555,
        416666,
        333333,
        277777,
    ],
    [
        10000000000,
        5000000000,
        100000000,
        33333333,
        50000000,
        25000000,
        16666666,
        12500000,
        10000000,
        8333333,
        5555555,
        4166666,
        3333333,
        2777777,
    ],
    [
        50000000000,
        1000000000,
        333333333,
        500000000,
        250000000,
        166666666,
        125000000,
        100000000,
        83333333,
        55555555,
        41666666,
        33333333,
        27777777,
    ],
]
# budgies = [ [None],
#
#             ]
big_budge = [
    [
        1000000000,
        333333333,
        250000000,
        166666666,
        125000000,
        83333333,
        55555555,
        41666666,
        27777777,
    ]
]
budge_abs, budge_res, budge_val = load_multiple_files_budget(
    [4],
    [35],
    30,
    nlsp_list=[1, 2, 3],
    N=10,
    Nts=100,
    max_copies=6,
    tags=["_nlsw1"],
    density_matrices=False,
    shots=budgies,
    folders=folders,
    train=True,
    budgets=[0, 5, 6, 7, 8, 9, 10, 11],
    train_use=100,
)
NLSP3NLSW1_abs = budge_abs
NLSP3NLSW0_abs, _, _ = load_multiple_files_budget(
    [4],
    [35],
    30,
    nlsp_list=[1, 2, 3],
    N=10,
    Nts=100,
    max_copies=6,
    tags=["_nlsw0"],
    density_matrices=False,
    shots=budgies,
    folders=folders,
    train=True,
    budgets=[0, 5, 6, 7, 8, 9, 10, 11],
    train_use=100,
)
# NLSP2FULL_abs  ,_,_ = load_multiple_files_budget([4,4,6,6,8,8,10],[4,4*16,6,6*16,8,8*16,10],30,nlsp_list = [1,2],N=10,Nts=100,max_copies=6,tags=tags,
#                                density_matrices = False,shots=big_budge,folders = folders,train=True,budgets=[1],train_use=100)
# NLSP2HALF_abs  ,_,_ = load_multiple_files_budget([4,4,6,6,8,8,10],[4,4*16,6,6*16,8,8*16,10],30,nlsp_list = [1,2],N=10,Nts=100,max_copies=6,tags=tags,
#                                density_matrices = False,shots=big_budge,folders = folders,train=True,budgets=[1],train_use=50)
# NLSP2QUART_abs ,_,_ = load_multiple_files_budget([4,4,6,6,8,8,10],[4,4*16,6,6*16,8,8*16,10],30,nlsp_list = [1,2],N=10,Nts=100,max_copies=6,tags=tags,
#                                density_matrices = False,shots=big_budge,folders = folders,train=True,budgets=[1],train_use=25)
# NLSP3QUART_abs ,_,_ = load_multiple_files_budget([4,4,6,6,8,8,10],[4,4*16,6,6*16,8,8*16,10],30,nlsp_list = [1,2,3],N=10,Nts=100,max_copies=6,tags=tags,
#                                density_matrices = False,shots=big_budge,folders = folders,train=True,budgets=[1],train_use=25)


# budge_abs.to_pickle('./BIG_CONGLOMERATE/'+"SWAP_NOISE_train_error.pkl")
# budge_res.to_pickle('./BIG_CONGLOMERATE/'+"SWAP_NOISE_train_rescaling.pkl")
# budge_val.to_pickle('./BIG_CONGLOMERATE/'+"SWAP_NOISE_train_full_results.pkl")
# NLSP2HALF_abs  = NLSP2HALF_abs.assign(description =  '2nlsp_half')
# NLSP2FULL_abs  = NLSP2FULL_abs.assign(description =     '2nlsp_full')
# NLSP2QUART_abs = NLSP2QUART_abs.assign(description = '2nlsp_quar')
NLSP3NLSW1_abs = NLSP3NLSW1_abs.assign(description="SWAP_NOISE(half)")
# NLSP3QUART_abs = NLSP3QUART_abs.assign(description = '3nlsp_quar')
NLSP3NLSW0_abs = NLSP3NLSW0_abs.assign(description="NORMAL(half)")

# NLSP2HALF_abs ,
# NLSP2FULL_abs ,
# NLSP2QUART_abs,
# NLSP3QUART_abs,
# ,
everything = pd.concat([NLSP3NLSW1_abs, NLSP3NLSW0_abs])
# %%
everything = pd.read_pickle(base_folder + "EVERYTHING_swap_noise.pkl")
# %%
everything_full = everything
everything_full.to_pickle(base_folder + "Swap_noise_complete_FULL_p35.pkl")
# %%
# budge_abs = pd.read_pickle('./BIG_CONGLOMERATE/'+"LINGROW2nlsp_abs_0filteredfull_error.pkl")
# everything_full = pd.read_pickle(base_folder+"EVERYTHING_FULL_more_Data_filt.pkl")

t_swap = filter_budget(
    everything, [0, 10**5, 10**6, 10**7, 10**8, 10**9, 10**10, 10**11]
)
t_swap["volume"] = t_swap.depth * t_swap.Qubits * t_swap.Qubits
t_swap["g"] = t_swap.depth // t_swap.Qubits
# %%

t_swap = pd.read_pickle(base_folder + "SWAP_DATA_v2_FULL_p35.pkl")

# t = pd.read_pickle('./BIG_CONGLOMERATE/'+"EVERYTHING_no_swap_nosie_filtered.pkl")
# %%
scaled_SWAP_Results = []
for i in [0]:
    for QUBIT, g in zip([4], [5]):
        temp = budge_abs.query(
            f"abs_error != 0 & nlsp == 1 & budget == {i} & Qubits == {QUBIT} & depth == {g}"
        ).copy()
        temp["scaled_error"] = (
            temp["abs_error"]
            / budge_abs.query(
                f'abs_error != 0 & nlsp == 1 & budget == {i} & copies == 1 & type == "VD_LinGrow" & Qubits == {QUBIT} & depth == {g}'
            )["abs_error"].mean()
        )
        scaled_SWAP_Results.append(temp)
scaled_swap_panda = pd.concat(scaled_SWAP_Results)

# %%
# %%
scaled_panda = pd.read_pickle(base_folder + "EVERYTHING_FULL_more_files_filtered.pkl")
# %%

t_swap["type"] = t_swap["type"].str.replace(r"_nlsw0", "")
t_swap["type"] = t_swap["type"].str.replace(r"_nlsw1", "")
ticker = t_swap.type.unique().tolist() + ["Noisy"]
colors = sns.color_palette("tab10", n_colors=len(ticker))  # get a number of colors
cmap = dict(zip(ticker, colors))  # zip values to colors
marker_map = {
    "VD": "s",
    "ZNE": "*",
    "CDR": "h",
    "vnCDR": "^",
    "vnCGVD": "p",
    "Noisy": "o",
    "CGVD": "v",
}
cmap = {
    "VD": "green",
    "ZNE": "red",
    "CDR": "pink",
    "vnCDR": "purple",
    "vnCGVD": "blue",
    "Noisy": "black",
    "CGVD": "orange",
}
line_map = {
    "VD": (4, 1.5),
    "ZNE": (4, 1.5),
    "CDR": (4, 1.5),
    "vnCDR": (4, 1.5),
    "vnCGVD": (4, 1.5),
    "Noisy": "",
    "CGVD": (4, 1.5),
}

# %%
tags = ["_nlsw0", "_nlsw1"]
title_dict = {"_nlsw0": "0 swap noise", "_nlsw1": "some swap noise"}
for tag in tags:
    for QUBIT, g in zip([4], [35]):
        fig, axes = plt.subplots(2, 3, figsize=(15, 10), sharey="row", sharex=True)
        axes = axes.flatten()
        dfs = []
        for i in [10**8, 10**9, 10**10, 10**11]:
            df1 = t_swap.query(
                f'budget=={i}    & Qubits == {QUBIT} & depth == {g} & abs_error != 0 & description.str.startswith("3") & description.str.endswith("half")'
            )
            df2 = t_swap.query(
                f'budget=={i*10}    & Qubits == {QUBIT} & depth == {g} & abs_error != 0 & ~(description.str.startswith("3") & description.str.endswith("half"))'
            )
            dfs.append(pd.concat([df1, df2]))
        sns.lineplot(
            data=t_swap.query(
                f'budget==0 &   Qubits == {QUBIT} & depth == {g} & abs_error != 0      &  type.str.endswith("{tag}")& type.str.startswith("vn")'
            ),
            ax=axes[0],
            x="copies",
            y="abs_error",
            style="type",
            hue="description",
            estimator="mean",
            markers=marker_map,
            err_style=None,
            palette=cmap,
        )
        sns.lineplot(
            data=t_swap.query(
                f'budget==0 &   Qubits == {QUBIT} & depth == {g} & abs_error != 0      &  type.str.endswith("{tag}")& type.str.startswith("vn")'
            ),
            ax=axes[3],
            x="copies",
            y="abs_error",
            style="type",
            hue="description",
            estimator="mean",
            markers=marker_map,
            err_style=None,
            palette=cmap,
        )

        sns.lineplot(
            data=dfs[0].query(
                f'   Qubits == {QUBIT} & depth == {g} & abs_error != 0      &  type.str.endswith("{tag}")& type.str.startswith("vn")'
            ),
            ax=axes[0 + 1],
            x="copies",
            y="abs_error",
            style="type",
            hue="description",
            estimator="mean",
            markers=marker_map,
            err_style=None,
            palette=cmap,
        )
        sns.lineplot(
            data=dfs[1].query(
                f'    Qubits == {QUBIT} & depth == {g} & abs_error != 0 &  type.str.endswith("{tag}") & type.str.startswith("vn")'
            ),
            ax=axes[1 + 1],
            x="copies",
            y="abs_error",
            style="type",
            hue="description",
            estimator="mean",
            markers=marker_map,
            err_style=None,
            palette=cmap,
        )
        sns.lineplot(
            data=dfs[2].query(
                f' Qubits == {QUBIT} & depth == {g} & abs_error != 0 &  type.str.endswith("{tag}") & type.str.startswith("vn")'
            ),
            ax=axes[2 + 2],
            x="copies",
            y="abs_error",
            style="type",
            hue="description",
            estimator="mean",
            markers=marker_map,
            err_style=None,
            palette=cmap,
        )
        sns.lineplot(
            data=dfs[3].query(
                f' Qubits == {QUBIT} & depth == {g} & abs_error != 0 &  type.str.endswith("{tag}") & type.str.startswith("vn")'
            ),
            ax=axes[3 + 2],
            x="copies",
            y="abs_error",
            style="type",
            hue="description",
            estimator="mean",
            markers=marker_map,
            err_style=None,
            palette=cmap,
        )
        for i in axes:
            i.set_yscale("log")
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
        fig.suptitle(
            f"different noise level and training data at {title_dict[tag]} Q={QUBIT} g={g}"
        )
        # plt.savefig(base_folder+f'/PLOTS/SWAP{tag}_budgetsQ{QUBIT}g{g}.png')
        # plt.close(fig)
# %%
o = t_swap.query(
    'res_type == "abs_error" & copies == 1 & nlsp == 1 & abs_error != 0  & budget != 0   & type.str.startswith("VD")'
).copy()
o["type"] = "Noisy"
a = (
    t_swap.query(
        'res_type == "abs_error" & copies == 1 & nlsp == 1 & abs_error !=  0  & budget != 0   & type.str.startswith("ZNE")'
    )
    .copy()
    .reset_index(drop=True)
)
b = t_swap.query(
    'res_type == "abs_error" & copies == 3 &  nlsp == 1 & abs_error != 0  & budget != 0  & type.str.startswith("VD")'
).copy()
c = t_swap.query(
    'res_type == "abs_error" & copies == 3 &  nlsp == 1 & abs_error != 0  & budget != 0  & type.str.startswith("vnCG")'
).copy()
d = t_swap.query(
    'res_type == "abs_error" & copies == 3 & nlsp == 1 & abs_error !=  0  & budget != 0   & type.str.startswith("vnCD")'
).copy()
e = t_swap.query(
    'res_type == "abs_error" & copies == 5 &  nlsp == 1 & abs_error != 0  & budget == 10**10 & budget != 0  & type.str.startswith("vnCG")'
).copy()
f = t_swap.query(
    'res_type == "abs_error" & copies == 5 & nlsp == 1 & abs_error !=  0  & budget == 10**9 & budget != 0   & type.str.startswith("vnCD")'
).copy()

d["shots"] = 2 * d["shots"]
c["shots"] = 2 * c["shots"]
c["budget"] = c["budget"] // 10
d["budget"] = d["budget"] // 10
c["shots"] = c["shots"] // 10
d["shots"] = d["shots"] // 10
e["budget"] = e["budget"] // 10
f["shots"] = f["shots"] // 10
e["shots"] = 2 * e["shots"]
f["shots"] = 2 * f["shots"]

selection = pd.concat([a, b, c, d, o, e, f])
# %%
o = t_swap.query(
    'res_type == "abs_error" & copies == 1 & nlsp == 1 & abs_error != 0  & budget != 0   & type.str.startswith("VD")'
).copy()
o["type"] = "Noisy"
a = (
    t_swap.query(
        'res_type == "abs_error" & copies == 1 & nlsp == 1 & abs_error !=  0  & budget != 0   & type.str.startswith("ZNE")'
    )
    .copy()
    .reset_index(drop=True)
)
b = t_swap.query(
    'res_type == "abs_error" & copies == 2 &  nlsp == 1 & abs_error != 0  & budget != 0  & type.str.startswith("VD")'
).copy()
c = t_swap.query(
    'res_type == "abs_error" & copies == 3 &  nlsp == 1 & abs_error != 0  & budget != 0  & type.str.startswith("vnCG")'
).copy()
d = t_swap.query(
    'res_type == "abs_error" & copies == 3 & nlsp == 1 & abs_error !=  0  & budget != 0   & type.str.startswith("vnCD")'
).copy()

c["budget"] = c["budget"]
d["budget"] = d["budget"]
c["shots"] = c["shots"]
d["shots"] = d["shots"]
e["budget"] = e["budget"]
f["shots"] = f["shots"]

selection = pd.concat([a, b, c, d, o])
# %%
o = t_swap.query(
    'res_type == "abs_error" & nlsp == 1 & abs_error != 0  & budget != 0   & type.str.startswith("VD")'
).copy()
o["type"] = "Noisy"
a = (
    t_swap.query(
        'res_type == "abs_error" & nlsp == 1 & abs_error !=  0  & budget != 0   & type.str.startswith("ZNE")'
    )
    .copy()
    .reset_index(drop=True)
)
b = t_swap.query(
    'res_type == "abs_error" &  nlsp == 1 & abs_error != 0  & budget != 0  & type.str.startswith("VD")'
).copy()
c = t_swap.query(
    'res_type == "abs_error" &  nlsp == 1 & abs_error != 0  & budget != 0  & type.str.startswith("vnCG")'
).copy()
d = t_swap.query(
    'res_type == "abs_error" & nlsp == 1 & abs_error !=  0  & budget != 0   & type.str.startswith("vnCD")'
).copy()
e = t_swap.query(
    'res_type == "abs_error" &  nlsp == 1 & abs_error != 0  & budget == 10**10 & budget != 0  & type.str.startswith("vnCG")'
).copy()
f = t_swap.query(
    'res_type == "abs_error" & nlsp == 1 & abs_error !=  0  & budget == 10**9 & budget != 0   & type.str.startswith("vnCD")'
).copy()

d["shots"] = 2 * d["shots"]
c["shots"] = 2 * c["shots"]
c["budget"] = c["budget"] // 10
d["budget"] = d["budget"] // 10
c["shots"] = c["shots"] // 10
d["shots"] = d["shots"] // 10
e["budget"] = e["budget"] // 10
f["shots"] = f["shots"] // 10
e["shots"] = 2 * e["shots"]
f["shots"] = 2 * f["shots"]

selection = pd.concat([a, b, c, d, o, e, f])
# %%
o = t_swap.query(
    'res_type == "abs_error" & nlsp == 1 & abs_error != 0  & budget != 0   & type.str.startswith("VD")'
).copy()
o["type"] = "Noisy"
a = (
    t_swap.query(
        'res_type == "abs_error" & nlsp == 1 & abs_error !=  0  & budget != 0   & type.str.startswith("ZNE")'
    )
    .copy()
    .reset_index(drop=True)
)
b = t_swap.query(
    'res_type == "abs_error" &  nlsp == 1 & abs_error != 0  & budget != 0  & type.str.startswith("VD")'
).copy()
c = t_swap.query(
    'res_type == "abs_error" &  nlsp == 1 & abs_error != 0  & budget != 0  & type.str.startswith("vnCG")'
).copy()
d = t_swap.query(
    'res_type == "abs_error" & nlsp == 1 & abs_error !=  0  & budget != 0   & type.str.startswith("vnCD")'
).copy()
e = t_swap.query(
    'res_type == "abs_error" &  nlsp == 1 & abs_error != 0  & budget == 10**10 & budget != 0  & type.str.startswith("vnCG")'
).copy()
f = t_swap.query(
    'res_type == "abs_error" & nlsp == 1 & abs_error !=  0  & budget == 10**9 & budget != 0   & type.str.startswith("vnCD")'
).copy()

d["shots"] = 2 * d["shots"]
c["shots"] = 2 * c["shots"]
c["budget"] = c["budget"]
d["budget"] = d["budget"]
c["shots"] = c["shots"]
d["shots"] = d["shots"]
e["budget"] = e["budget"]
f["shots"] = f["shots"]
e["shots"] = 2 * e["shots"]
f["shots"] = 2 * f["shots"]

selection = pd.concat([a, b, c, d, o, e, f])


# %%
o = t_swap.query(
    ' nlsp == 1 & abs_error != 0   & description == "3nlsp_full"  & type.str.startswith("VD")'
).copy()
o["type"] = "Noisy"
a = (
    t_swap.query(
        ' nlsp == 1 & abs_error != 0   & description == "3nlsp_full"  & type.str.startswith("ZNE")'
    )
    .copy()
    .reset_index(drop=True)
)
b = t_swap.query(
    '  nlsp == 1 & abs_error != 0   & description == "3nlsp_full"  & type.str.startswith("VD")'
).copy()
c = t_swap.query(
    '  nlsp == 1 & abs_error != 0   & description == "3nlsp_half"  & type.str.startswith("vnCG")'
).copy()
d = t_swap.query(
    ' nlsp == 1 & abs_error != 0   & description == "3nlsp_half"  & type.str.startswith("vnCD")'
).copy()
i = t_swap.query(
    '  nlsp == 1 & abs_error != 0   & description == "3nlsp_half"  & type.str.startswith("CG")'
).copy()
j = t_swap.query(
    'copies==1 & nlsp == 1 & abs_error != 0   & description == "3nlsp_half"  & type.str.startswith("CD")'
).copy()
d["shots"] = 2 * d["shots"]
c["shots"] = 2 * c["shots"]
i["shots"] = 2 * i["shots"]
j["shots"] = 2 * j["shots"]
e = t_swap.query(
    '  nlsp == 1 & abs_error != 0   & description == "3nlsp_full"  & type.str.startswith("vnCG")'
).copy()
f = t_swap.query(
    ' nlsp == 1 & abs_error != 0   & description == "3nlsp_full"  & type.str.startswith("vnCD")'
).copy()
g = t_swap.query(
    '  nlsp == 1 & abs_error != 0   & description == "3nlsp_full"  & type.str.startswith("CG")'
).copy()
h = t_swap.query(
    'copies==1 & nlsp == 1 & abs_error != 0   & description == "3nlsp_full"  & type.str.startswith("CD")'
).copy()

idnex_mis = a.query("copies == 1  & budget==10**9 & Qubits == 6 & g == 1").index
idnex_rep = a.query("copies == 1  & budget==10**8 & Qubits == 6 & g == 1").index
a.loc[idnex_mis, "abs_error"] = list(a.iloc[idnex_rep]["abs_error"])
idnex_mis = a.query("copies == 1  & budget==10**9 & Qubits == 10 & g == 1").index
idnex_rep = a.query("copies == 1  & budget==10**8 & Qubits == 10 & g == 1").index
a.loc[idnex_mis, "abs_error"] = list(a.iloc[idnex_rep]["abs_error"])

c["budget"] = c["budget"] // 10
d["budget"] = d["budget"] // 10
c["shots"] = c["shots"] // 10
d["shots"] = d["shots"] // 10
j["budget"] = j["budget"] // 10
i["budget"] = i["budget"] // 10
i["shots"] = i["shots"] // 10
j["shots"] = j["shots"] // 10

selection = pd.concat([a, b, c, d, o, e, f, g, h, i, j])
selection["type"] = selection["type"].str.replace(r"_LinGrow", "")
# %%
#### with inf
o = t_swap.query(
    'res_type == "abs_error" & copies == 1 & nlsp == 1 & abs_error != 0 & budget !=0   & type.str.startswith("VD")'
).copy()
o["type"] = "Noisy"
a = (
    t_swap.query(
        'res_type == "abs_error" & copies == 1 & nlsp == 1 & abs_error !=  0  & budget !=0    & type.str.startswith("ZNE")'
    )
    .copy()
    .reset_index(drop=True)
)
b = t_swap.query(
    'res_type == "abs_error" & copies == 3 &  nlsp == 1 & abs_error != 0  & budget !=0    & type.str.startswith("VD")'
).copy()
c = t_swap.query(
    'res_type == "abs_error" & copies == 3 &  nlsp == 1 & abs_error != 0  & budget !=0   & type.str.startswith("vnCG")'
).copy()
d = t_swap.query(
    'res_type == "abs_error" & copies == 3 & nlsp == 1 & abs_error !=  0  & budget !=0     & type.str.startswith("vnCD")'
).copy()
e = t_swap.query(
    'res_type == "abs_error" & copies == 5 &  nlsp == 1 & abs_error != 0  & budget == 10**10 & budget != 0  & type.str.startswith("vnCG")'
).copy()
f = t_swap.query(
    'res_type == "abs_error" & copies == 5 & nlsp == 1 & abs_error !=  0  & budget == 10**9 & budget != 0   & type.str.startswith("vnCD")'
).copy()

o_ = t_swap.query(
    'res_type == "abs_error" & copies == 1 & nlsp == 1 & abs_error != 0 & budget ==0   & type.str.startswith("VD")'
).copy()
o_["type"] = "Noisy"
a_ = (
    t_swap.query(
        'res_type == "abs_error" & copies == 1 & nlsp == 1 & abs_error !=  0  & budget ==0    & type.str.startswith("ZNE")'
    )
    .copy()
    .reset_index(drop=True)
)
b_ = t_swap.query(
    'res_type == "abs_error" & copies == 3 &  nlsp == 1 & abs_error != 0  & budget ==0    & type.str.startswith("VD")'
).copy()
c_ = t_swap.query(
    'res_type == "abs_error" & copies == 4 &  nlsp == 1 & abs_error != 0  & budget ==0   & type.str.startswith("vnCG")'
).copy()
d_ = t_swap.query(
    'res_type == "abs_error" & copies == 3 & nlsp == 1 & abs_error !=  0  & budget ==0     & type.str.startswith("vnCD")'
).copy()
o_["budget"] = 10**12
a_["budget"] = 10**12
b_["budget"] = 10**12
c_["budget"] = 10**12
d_["budget"] = 10**12

d["shots"] = 2 * d["shots"]
c["shots"] = 2 * c["shots"]
c["budget"] = c["budget"]
d["budget"] = d["budget"]
c["shots"] = c["shots"]
d["shots"] = d["shots"]
e["budget"] = e["budget"]
f["shots"] = f["shots"]
e["shots"] = 2 * e["shots"]
f["shots"] = 2 * f["shots"]

selection = pd.concat([a, b, c, d, o, e, f, a_, b_, c_, d_, o_])

# %%

g = sns.FacetGrid(selection, col="description", sharey=True)
g.map_dataframe(
    sns.lineplot,
    x="budget",
    y="abs_error",
    hue="type",
    style="type",
    dashes=line_map,
    markers=marker_map,
    err_style=None,
    palette=cmap,
).set(
    yscale="log",
    xscale="log",
    xlim=[10**5, 10**12 + 10**8],
    xticks=[10**5, 10**6, 10**7, 10**8, 10**9, 10**10, 10**11, 10**12],
    xticklabels=[
        r"$10^6$",
        r"$10^6$",
        r"$10^7$",
        r"$10^8$",
        r"$10^9$",
        r"$10^{10}$",
        r"$10^{11}$",
        r"$\infty$",
    ],
)
g.add_legend()
g.set_axis_labels("Budget", r"Mean $\langle \sigma^1_z \rangle$ absolute error")
g.tight_layout()
# g.fig.delaxes(g.axes[-1, -1])
# g.axes[0][0].set_ylabel('')
# g.axes[1][0].yaxis.set_label_coords(-0.125,1.08)
# for ax in g.axes[0]:
#     ax.set_ylim([10**-4,10**-1])
# for ax in g.axes[1]:
#     ax.set_ylim([10**-2,10**0])
# plt.savefig(base_folder+'abs_error_per_budget.pdf')
# %%

o = t_swap.query(
    'nlsp == 1 & abs_error != 0  & budget > 10**9 & description == "3nlsp_full" & copies == 1 & type.str.startswith("VD")'
).copy()
o["type"] = "Noisy"
a = t_swap.query(
    'nlsp == 1 & abs_error != 0  & budget > 10**9 & description == "3nlsp_full" & copies == 1 & type.str.startswith("ZNE")'
).copy()
b = t_swap.query(
    'nlsp == 1 & abs_error != 0  & budget > 10**9 & description == "3nlsp_full" & copies == 3 & type.str.startswith("VD")'
).copy()
c = t_swap.query(
    'nlsp == 1 & abs_error != 0  & budget > 10**10 & description == "3nlsp_half" & copies == 3 & type.str.startswith("vnCG")'
).copy()
d = t_swap.query(
    'nlsp == 1 & abs_error != 0  & budget > 10**10 & description == "3nlsp_half" & copies == 1 & type.str.startswith("vnCD")'
).copy()

c["budget"] = c["budget"] // 10
d["budget"] = d["budget"] // 10

selection = pd.concat([a, b, c, d, o])

g = sns.FacetGrid(selection, col="Qubits", row="g")
g.map_dataframe(sns.boxplot, x="budget", y="abs_error", hue="type", palette=cmap).set(
    yscale="log"
)
g.add_legend(bbox_to_anchor=(0.85, 0.3))
g.set_axis_labels("Budget", r"Mean $\langle \sigma^1_z \rangle$ absolute error")
# g.fig.delaxes(g.axes[-1, -1])
# g.axes[0][0].set_ylabel('')
# g.axes[1][0].yaxis.set_label_coords(-0.25,1.08)
plt.savefig(base_folder + "abs_error_over_budget.pdf")
# %%
import csv

with open("dataframe_avg.txt", "w") as csv_file:
    csv_writer = csv.writer(csv_file, delimiter=",")
    line_count = 0
    csv_writer.writerow(
        [
            "type",
            "budget",
            "qubit",
            "g",
            "value",
            "shots",
            "copies",
            "description",
            "result_type",
        ]
    )
    for thing in selection.type.unique():
        print(thing)
        for nlsp in [1]:
            for budget in [10**5, 10**6, 10**7, 10**8, 10**9, 10**10]:
                for q, g in zip([4, 6, 8, 10, 4, 6, 8], [1, 1, 1, 1, 16, 16, 16]):
                    for result_type in ["abs_error", "mit_result", "exact_value"]:
                        x = selection.query(
                            f" res_type == '{result_type}' &budget == {budget} & type == '{thing}' & nlsp == {nlsp} & Qubits == {q} & g == {g}"
                        )["abs_error"].mean()
                        frame = selection.query(
                            f" budget == {budget} & type == '{thing}' & nlsp == {nlsp} & Qubits == {q} & g == {g}"
                        )
                        try:
                            desc = str(frame["description"].iloc[0])
                            sho = str(frame["shots"].iloc[0])
                            cop = str(frame["copies"].iloc[0])
                        except:
                            desc = "no"
                            sho = "no"
                            cop = "no"
                        csv_writer.writerow(
                            [thing, budget, q, g, x, sho, cop, desc, result_type]
                        )

description = "3nlsp_half"
title_dict = {
    "_LinGrow": "Linear Growth"
}  # ,'_RNM_e1':'Scaled Noise (e_r=1)','_RNM_e0p1':'Scaled Noise (e_r=0.1)'}
for budget in [10**7, 10**8, 10**9, 10**10]:
    fig, axes = plt.subplots(2, 4, figsize=(15, 10), sharey=False, sharex="col")
    #    axes=axes.flatten():q
    plot_data = []
    plot_data.append(
        scaled_panda.query(
            f'(description == "{description}" & budget=={budget}|budget==0)  & copies==1  & nlsp == 1  & abs_error != 0  &  (type.str.startswith("ZNE")|type.str.startswith("CDR")|type.str.startswith("vnCD"))'
        )
    )
    plot_data.append(
        scaled_panda.query(
            f'(description == "{description}" & budget=={budget}|budget==0)  & copies==3  & nlsp == 1  & abs_error != 0  &  (type.str.startswith("VD")|type.str.startswith("vnCG")|type.str.startswith("CG"))'
        )
    )

    plot_data = pd.concat(plot_data)
    for row, tag in enumerate(["_LinGrow"]):
        for column, g in enumerate([1, 16]):
            sns.lineplot(
                data=plot_data.query(
                    f'    budget=={budget} & nlsp == 1 &  abs_error != 0  &  type.str.endswith("{tag}") & g == {g}'
                ),
                ax=axes[column, row],
                x="Qubits",
                y="scaled_error",
                hue="type",
                style="type",
                estimator="mean",
                markers=marker_map,
                err_style=None,
                palette=cmap,
                color={"VD_LinGrow": "pink", "vnCDR_LinGrow": "#742802"},
            )
            axes[column, row].set_yscale("log")
            axes[column, row].legend([], [], frameon=False)
            axes[column, row].set_title(f"Model: {title_dict[tag]}, g = {g}")

    handles, labels = axes[0, 0].get_legend_handles_labels()
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
    if budget != 0:
        fig.suptitle(
            f"More shots CG/vnPerformance over qubits and depth at budget 10^{int(np.log10(budget))}  "
        )
        plt.savefig(f"./PLOTS/nlsp3_full-10^{int(np.log10(budget))}screens.png")

    else:
        fig.suptitle("Performance over qubits and depth at infinite shot limit.")
        plt.savefig("./PLOTS/nlsp3_full-budget_inf_shots_screens.png")
    plt.close(fig)


# %%
base_folder = "D:/Databases/VD_CDR/"
fig, axes = plt.subplots(6, 2, figsize=(15, 15), sharey=False, sharex="col")
#    axes=axes.flatten():q
description = "hlaf_v_normal"
t_swap["type"] = t_swap["type"].str.replace(r"_LinGrow", "")

# plot_data=[]
# plot_data.append(scaled_panda.query(f'description == "{description}" & copies==1  & nlsp == 1  & abs_error != 0  &  type.str.startswith("ZNE")'))
# plot_data.append(scaled_panda.query(f'description == "{description}" &  copies==3  & nlsp == 1  & abs_error != 0  &  (type.str.startswith("VD")|type.str.startswith("vn")|type.str.startswith("CG"))'))
# plot_data=pd.concat(plot_data)

for row, budg in enumerate([10**5, 10**6, 10**7, 10**8, 10**9, 10**10]):
    for column, g in enumerate([1, 16]):
        plot_data = []
        plot_data.append(
            t_swap.query(
                f'budget == {budg} & description == "3nlsp_full" & copies==1  & nlsp == 1  & abs_error != 0  &  (type.str.startswith("ZNE"))'
            )
        )
        plot_data.append(
            t_swap.query(
                f'budget == {budg} & description == "3nlsp_full" & copies==3  & nlsp == 1  & abs_error != 0  &  (type.str.startswith("VD"))'
            )
        )
        plot_data.append(
            t_swap.query(
                f'budget == {budg*10} & description == "3nlsp_half" &  copies==3  & nlsp == 1  & abs_error != 0  &  type.str.startswith("vnCG")'
            )
        )
        plot_data.append(
            t_swap.query(
                f'budget == {budg*10} & description == "3nlsp_half" &  copies==1  & nlsp == 1  & abs_error != 0  &  type.str.startswith("vnCD")'
            )
        )
        o = t_swap.query(
            'abs_error != 0  & budget != 0 & description == "3nlsp_full" & copies == 1 & type.str.startswith("VD")'
        ).copy()
        o["type"] = "Noisy"
        plot_data.append(o)
        plot_data = pd.concat(plot_data)

        sns.lineplot(
            data=plot_data.query(f" g == {g}"),
            ax=axes[row, column],
            x="Qubits",
            y="abs_error",
            hue="type",
            style="type",
            estimator="mean",
            markers=marker_map,
            err_style=None,
            palette=cmap,
            dashes=line_map,
        )
        axes[row, column].set_yscale("log")
        axes[row, column].legend([], [], frameon=False)
        # axes[row,column].set_ylim([10**-2,10**1])
        # axes[row,1].set_yticklabels([])
        axes[row, 1].set_ylabel("")
        axes[row, 0].set_ylabel(r"$\overline{\langle \sigma_z \rangle}$ abs error")

        try:
            axes[row, column].set_title(f"Budget: 10^{int(np.log10(budg))}, g = {g}")
        except:
            axes[row, column].set_title(f"Budget: Infinite, g = {g}")


handles, labels = axes[-2][-2].get_legend_handles_labels()
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
if budg != 0:
    fig.suptitle(f" Comparison between budgets at different techniques")
    plt.savefig(base_folder + f"/PLOTS/A{description}_lineplot_budget_screens.png")

else:
    fig.suptitle("Comparison between budgets")
    plt.savefig(base_folder + f"/PLOTS/A{description}_lineplot_budget_screens.png")
plt.close(fig)

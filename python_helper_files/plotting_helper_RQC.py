import pandas as pd
import seaborn as sns
from pathlib import Path
from post_processing_multi import *

tags = [""]
# predefined budgets
budget_shots = [
    [None],
    {
        50000,
        1000,
        333,
        500,
        250,
        166,
        125,
        100,
        83,
        55,
        41,
        33,
        27,
        500000,
        10000,
        3333,
        5000,
        2500,
        1666,
        1250,
        1000,
        833,
        555,
        416,
        333,
        277,
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
    },
]


def relative_shots(training_sets, noise_levels, max_copies, shots_budget):
    """Returns relative shots with less shots for COI when shots budget is high enough.

    Args:
        training_sets (int): number of training sets.
        noise_levels (int):  number of  noise levels
        max_copies (int): number of copies considred
        shots_budget (int): total budget

    Returns:
        pandas DF: dataframe with shot budgets
    """
    if shots_budget != 0:
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
                "UNITED": [
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
                "UNITED": [None] * 6,
                "copies": list(range(1, max_copies + 1)),
            }
        )
    return pd.melt(storage, id_vars=("copies"))


def relative_shots_table(
    training_sets, noise_levels, max_copies, shots_budget, swap_levels
):
    """Returns relative shots used by techniques with no optimizations or considerations. Includes
    swap nosie.

    Args:
        training_sets (int): number of training sets.
        noise_levels (int):  number of  noise levels
        max_copies (int): number of copies considred
        shots_budget (int): total budget
        swap_levels (int): number of swap noise levels

    Returns:
        pandas DF: dataframe with shot budgets
        ratio: adjusted budget from full
    """
    ratio = int(shots_budget // int((100 // training_sets) * (max_copies - 1)))
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
                "UNITED": [vnCDR_shots]
                + [
                    x // noise_levels // 2
                    for x in [
                        base_cost // (copy + 2)
                        for copy, base_cost in enumerate(
                            [CGVD_shots] * (max_copies - 1)
                        )
                    ]
                ],
                "UNITED+": [vnCDR_shots // swap_levels]
                + [
                    x // noise_levels // 2 // swap_levels
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
                "UNITED": None,
                "copies": list(range(1, max_copies + 1)),
            }
        )
    return storage, ratio


def single_qubit_gate_length(layers, qubits):
    return layers * 6 * (qubits - 1)


def two_qubit_gate_length(layers, qubits):
    return layers * (qubits - 1)


def circuit_length(layers, qubits):
    """calculates circuit length

    Args:
        layers (int): number of layers in RQC
        qubits (int): number of qubtis in system

    Returns:
        int: decomposed circuit depth
    """
    e1 = layers * 6 * (qubits - 1)
    e2 = layers * (qubits - 1)
    return e1 + e2


def total_error_rates_normal(layers, qubits, factor=1):
    """Generates total error rate in standard numerical experiment.

    Args:
        layers (int): number of layers of RQC.
        qubits (int): number of qubits
        factor (float, optional): noise scaling factor. Defaults to 1.

    Returns:
        float: total error rate.
    """
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
    """Loads multiple files and processes all the data through the various techniques.

    Args:
        Q (list): list of qubits
        p (list): depth list, same length as Q
        num_seeds (int): number of seeds considered
        nlsp_list (list, optional): list of noise levels. Defaults to [1,2,3].
        N (int, optional): non clifford gates in simulation. Defaults to 10.
        Nts (int, optional): training set size. Defaults to 100.
        max_copies (int, optional): max number of copies in VD. Defaults to 6.
        tags (list, optional): list of additional tags. Defaults to [''].
        density_matrices (bool, optional): Do you have density matrices?. Defaults to False.
        shots (list, optional): shots of finite shot simulations. Defaults to [None].
        folder (str, optional): folder to be analyzed. Defaults to ''.
        train (bool, optional): existance of training data. Defaults to True.

    Returns:
        (absolute_error_df, error_rescaling_df, mitigated_values_df) : Various dfs of processed and organized data.
    """
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
    """Loads multiple files and processes all the data through the various techniques, assigns budget as well.
    But this needs to be further filtered.

    Args:
        Q (list): list of qubits
        p (list): depth list, same length as Q
        num_seeds (int): number of seeds considered
        nlsp_list (list, optional): list of noise levels. Defaults to [1,2,3].
        N (int, optional): non clifford gates in simulation. Defaults to 10.
        Nts (int, optional): training set size. Defaults to 100.
        max_copies (int, optional): max number of copies in VD. Defaults to 6.
        tags (list, optional): list of additional tags. Defaults to [''].
        density_matrices (bool, optional): Do you have density matrices?. Defaults to False.
        shots (list[lists], optional): list of lists of shots of finite shot simulations. Defaults to [[None]].
        folder (str, optional): folder to be analyzed. Defaults to ''.
        train (bool, optional): existance of training data. Defaults to True.
        budgets (list): list with same number of entries as shots lists corresponding to relative budget.
        train_use (int): how many trainig sets do you want to use?

    Returns:
        (absolute_error_df, error_rescaling_df, mitigated_values_df) : Various dfs of processed and organized data.
    """

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
    dfs_absolute_error = pd.concat(dfs_absolute_error)
    dfs_error_rescaling = pd.concat(dfs_error_rescaling)
    dfs_standard_values = pd.concat(dfs_standard_values)
    return dfs_absolute_error, dfs_error_rescaling, dfs_standard_values


def filter_budget(df, budgets, training=100, max_copies=6, noise_levels=3):
    """Filters appropiate budget and shots for the technique in question.

    Args:
        df (result_df from load_multiple_files): Appropiate dataframe.
        budgets (list): list of budgets.
        training_use (int): how many training samples are actually used.
        noise_levels (int): maximum number of noise levels considered.
        copies (int): maximum number of copies considered.

    Returns:
        filtered_dataframe: dataframe with correct shot count data for appropiate technique
    """

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
                if (technique == "UNITED" or technique == "vnCDR") and training == 50:
                    df_temp_abs = df.query(
                        f'shots == {allowance} & type.str.startswith("{technique}") & copies == {copy_no}'
                    ).assign(budget=compute_ratio(budget))
                else:
                    df_temp_abs = df.query(
                        f'shots == {allowance} & type.str.startswith("{technique}") & copies == {copy_no}'
                    ).assign(budget=budget)
            filtered_dfs.append(df_temp_abs)
    filtered_dfs = pd.concat(filtered_dfs)
    return filtered_dfs


def plot_over_budget(df):
    zero_copy_methods = df.query(
        'abs_error > 0  & copies == 1 & nlsp==1 & description == "3nlsp_full" & res_type=="abs_error" & ( type == "ZNE" | type == "vnCDR")'
    )
    noisy = df.query(
        'abs_error > 0  & copies == 1 & nlsp==1 & description == "3nlsp_full" & res_type=="abs_error" & ( type=="VD")'
    )
    few_copy_methods = df.query(
        'abs_error >  0  & nlsp==1 & copies==2 & description == "3nlsp_full" & res_type=="abs_error" & ( type=="VD")'
    )
    many_copy_methods = df.query(
        'abs_error > 0  &nlsp==1  & copies==2& description == "3nlsp_full" & res_type=="abs_error" & ( type=="UNITED")'
    )
    noisy["type"] = "noisy"
    plot_df = pd.concat(
        [noisy, zero_copy_methods, few_copy_methods, many_copy_methods],
        axis=0,
        ignore_index=True,
    )
    fig = sns.relplot(
        data=plot_df.reset_index(),
        kind="line",
        x="budget",
        y="abs_error",
        hue="type",
        col="Qubits",
        col_wrap=2,
        style="copies",
        estimator="mean",
        markers=True,
        ci=None,
    ).set(yscale="log", xscale="log")
    ax = fig.axes[0]
    xticks = ax.get_xticks().tolist()
    for i in range(len(xticks) - 1):
        xticks[i] = rf"$10^{{{i+3}}}$"
    xticks[-3] = r"$\infty$"
    xticks[-4] = r"$\cdots$"
    ax.set_xticklabels(xticks)
    return fig


def compute_ratio(budget):
    table, ratio = relative_shots_table(50, 3, 6, budget, 1)
    return ratio


def load_rqc_data(
    qubit_list=[4, 4, 6, 6, 8, 8, 10], depths=[4, 4 * 16, 6, 6 * 16, 8, 8 * 16, 10]
):
    if Path("./RQC_runs/checkpoint_full_training.pq").is_file():
        everything = pd.concat([pd.read_parquet("./RQC_runs/checkpoint_full_training.pq"),pd.read_parquet("./RQC_runs/checkpoint_half_training.pq")])
    else:
        base_folder = "./RQC_runs/all_qubits/"
        folders = [base_folder]
        tags = ["_LinGrow"]
        budgies = [
            [None],
            {
                50000,
                1000,
                333,
                500,
                250,
                166,
                125,
                100,
                83,
                55,
                41,
                33,
                27,
                500000,
                10000,
                3333,
                5000,
                2500,
                1666,
                1250,
                1000,
                833,
                555,
                416,
                333,
                277,
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
            },
        ]
        qubit_list = [4, 4, 6, 6, 8, 8, 10]
        depths = [4, 4 * 16, 6, 6 * 16, 8, 8 * 16, 10]
        budge_abs, budge_res, budge_val = load_multiple_files_budget(
            qubit_list,
            depths,
            30,
            nlsp_list=[1, 2, 3],
            N=10,
            Nts=100,
            max_copies=6,
            tags=tags,
            density_matrices=False,
            shots=budgies,
            folders=folders,
            train=True,
            budgets=[0, 5],
            train_use=100,
        )
        NLSP3FULL_abs = budge_abs
        NLSP3HALF_abs, _, _ = load_multiple_files_budget(
            qubit_list,
            depths,
            30,
            nlsp_list=[1, 2, 3],
            N=10,
            Nts=100,
            max_copies=6,
            tags=tags,
            density_matrices=False,
            shots=budgies,
            folders=folders,
            train=True,
            budgets=[0, 5],
            train_use=50,
        )
        NLSP3FULL_abs = NLSP3FULL_abs.assign(description="3nlsp_full")
        NLSP3HALF_abs = NLSP3HALF_abs.assign(description="3nlsp_half")
        everything = pd.concat([NLSP3FULL_abs, NLSP3HALF_abs])
        everything["type"] = everything["type"].str.replace("_LinGrow", "")
        everything.query('description=="3nlsp_full"').to_parquet("./RQC_runs/checkpoint_full_training.pq")
        everything.query('description=="3nlsp_half"').to_parquet("./RQC_runs/checkpoint_half_training.pq")
    return everything


def filter_rqc_data(
    df, budgets=[10 ** 5, 10 ** 6, 10 ** 7, 10 ** 8, 10 ** 9, 10 ** 10, 10 ** 11, 0]
):
    if Path("./RQC_runs/filtered_data.pq").is_file():
        filtered_rqc = pd.read_parquet("./RQC_runs/filtered_data.pq")
    else:
        rqc_half_df = filter_budget(df.query('description=="3nlsp_half"'), budgets, 50)
        rqc_half_df["volume"] = (
            rqc_half_df.depth * rqc_half_df.Qubits * rqc_half_df.Qubits
        )
        rqc_half_df["g"] = rqc_half_df.depth // rqc_half_df.Qubits

        rqc_full_df = filter_budget(df.query('description=="3nlsp_full"'), budgets, 100)
        rqc_full_df["volume"] = (
            rqc_full_df.depth * rqc_full_df.Qubits * rqc_full_df.Qubits
        )
        rqc_full_df["g"] = rqc_full_df.depth // rqc_full_df.Qubits
        filtered_rqc = pd.concat([rqc_full_df, rqc_half_df])
        filtered_rqc.to_parquet("./RQC_runs/filtered_data.pq")
    return filtered_rqc


def figure_3(df, statistic_to_plot="mean", qubit=4, training="half"):
    """Plots figure 3, which is the performance of the techniques over budget. This is extended in the sense that it can return
    all qubit results in a panel."""
    df = df.query(f'Qubits=={qubit} & description == "3nlsp_{training}"')

    zero_copy_methods = df.query(
        'abs_error > 0  & copies == 1 & nlsp==1  & res_type=="abs_error" & ( type == "ZNE" )'
    )
    vnCDR = df.query(
        'type == "vnCDR"&abs_error > 0  & copies == 1 & nlsp==1  & res_type=="abs_error"'
    )
    noisy = df.query(
        'abs_error > 0  & copies == 1 & nlsp==1 & res_type=="abs_error" & ( type=="VD")'
    )
    few_copy_methods = df.query(
        'abs_error >  0  & nlsp==1 & copies==2 & res_type=="abs_error" & ( type=="VD")'
    )
    many_copy_methods = df.query(
        'abs_error > 0  &nlsp==1  & copies==3 & res_type=="abs_error" & ( type=="UNITED")'
    )

    noisy["type"] = "noisy"
    plot_df = pd.concat(
        [noisy, zero_copy_methods, few_copy_methods, many_copy_methods, vnCDR],
        axis=0,
        ignore_index=True,
    )
    fig = sns.relplot(
        data=plot_df.reset_index().query("budget>10**4&budget<10**11"),
        kind="line",
        x="budget",
        y="abs_error",
        hue="type",
        col="Qubits",
        row="g",
        estimator="mean",
        markers=True,
        ci=None,
    ).set(xscale="log", yscale="log")


def figure_2(df, statistic_to_plot="mean", like_paper=True, g=1):
    """"Plot figure 2 from the paper, which is the absoulte error over qubits at different budgets."""
    if like_paper:
        df = df.query(
            f"budget>10**4&budget<10**11&budget!=10**9&budget!=10**7 & g=={g}"
        )
    else:  # Basically will show all computed budgets
        df = df.query(f"budget>10**4&budget<10**11&g=={g}")
    zero_copy_methods = df.query(
        'abs_error > 0  & copies == 1 & nlsp==1 & description == "3nlsp_half" & res_type=="abs_error" & ( type == "ZNE" | type == "vnCDR")'
    )
    noisy = df.query(
        'abs_error > 0  & copies == 1 & nlsp==1 & description == "3nlsp_full" & res_type=="abs_error" & ( type=="VD")'
    )
    few_copy_methods = df.query(
        'abs_error >  0  & nlsp==1 & copies==2 & description == "3nlsp_full" & res_type=="abs_error" & ( type=="VD")'
    )
    many_copy_methods = df.query(
        'abs_error > 0  &nlsp==1  & copies==3& description == "3nlsp_half" & res_type=="abs_error" & ( type=="UNITED")'
    )
    noisy["type"] = "noisy"
    plot_df = pd.concat(
        [noisy, zero_copy_methods, few_copy_methods, many_copy_methods],
        axis=0,
        ignore_index=True,
    )
    fig = sns.relplot(
        data=plot_df.reset_index(),
        kind="line",
        x="Qubits",
        col="budget",
        col_wrap=2,
        y="abs_error",
        hue="type",
        estimator=statistic_to_plot,
        markers=True,
        ci=None,
    )
    fig.set(yscale="log")
    fig.set_ylabels(statistic_to_plot + " of absolute error")
    return fig


def load_raw_rqc_data():
    """Returns a dictionary of the coi data and training data with indices 'coi' and 'train'."""
    import os

    dir_path = "RQC_runs/all_qubits/"

    def quick_load(fi, dir_path):
        res = []
        files = []
        for file in os.listdir(dir_path):
            if file.endswith(".pq") and file.startswith(f"pandas_{fi}"):
                res.append(
                    pd.read_parquet(dir_path + file)
                    .assign(file=file)
                    .assign(data=fi)
                    .reset_index()
                )
                files.append(file)
        res = pd.concat(res)
        res = res.reset_index()
        res = res.reset_index()
        res = res.drop(["result_type"], axis=1)
        res["exact"] = pd.to_numeric(res["exact"])
        res["expectation"] = pd.to_numeric(res["expectation"])
        res["exact_abs"] = np.abs(res["exact"])
        return res

    train = quick_load("train", dir_path)
    coi = quick_load("COI", dir_path)
    return {"coi": coi, "train": train}


def figure_7(df, budget=10 ** 10, statistic_to_plot="mean", g=1):
    """Plots figure 7 for any budget and g value (depth scaling)"""
    df = df.query(f"budget=={budget} & g=={g}")
    few_copy_methods = df.query(
        'abs_error >  0  & nlsp==1  & description == "3nlsp_half" & res_type=="abs_error" & ( type=="VD")'
    )
    many_copy_methods = df.query(
        'abs_error > 0  &nlsp==1  & description == "3nlsp_half" & res_type=="abs_error" & ( type=="UNITED")'
    )
    plot_df = pd.concat([few_copy_methods], axis=0, ignore_index=True,)
    fig = sns.relplot(
        data=plot_df.reset_index(),
        kind="line",
        x="copies",
        col="g",
        y="abs_error",
        hue="Qubits",
        style="type",
        estimator=statistic_to_plot,
        ci=None,
    )
    fig.set(yscale="log", title=f"budget {budget}")
    fig.set_ylabels(statistic_to_plot + " of absolute error")
    return fig


def figure_8(coi, train, qubit=4, depth=4, g=None):
    """Plots the disctibution of training data and circuit of interest for desired qubit"""
    if g is not None:
        if g == 1 or g == 16:
            depth = int(qubit * g)
        else:
            print("g value nonexistent")
    npcoi = coi.query(f"file.str.contains('{qubit}p{depth}')")
    nptrain = train.query(f"file.str.contains('{qubit}p{depth}')")
    sns.displot(
        kind="hist",
        data=pd.concat([npcoi, nptrain]),
        x="exact",
        col="file",
        col_wrap=4,
        bins=25,
        facet_kws=dict(sharey=True),
    ).set(yscale="log")


def effect_of_training_set(df):
    zero_copy_methods = df.query(
        'abs_error > 0  & copies == 1 & nlsp==1  & res_type=="abs_error" & ( type == "ZNE" )'
    )
    vnCDR = df.query(
        'type == "vnCDR"&abs_error > 0  & copies == 1 & nlsp==1  & res_type=="abs_error"'
    )
    noisy = df.query(
        'abs_error > 0  & copies == 1 & nlsp==1 & res_type=="abs_error" & ( type=="VD")'
    )
    few_copy_methods = df.query(
        'abs_error >  0  & nlsp==1 & copies==2 & res_type=="abs_error" & ( type=="VD")'
    )
    many_copy_methods = df.query(
        'abs_error > 0  &nlsp==1  & copies==3 & res_type=="abs_error" & ( type=="UNITED")'
    )
    plot_df = pd.concat([many_copy_methods, vnCDR], axis=0, ignore_index=True,)
    fig = sns.relplot(
        data=plot_df.reset_index().query("(budget==0|budget==10**10)&g==1"),
        kind="line",
        y="abs_error",
        x="Qubits",
        hue="type",
        style="description",
        col="budget",
        estimator="max",
        markers=True,
        ci=None,
    )
    fig.set(yscale="log")
    return fig

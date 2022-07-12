import pandas as pd
import seaborn as sns
from post_processing_multi import *

tags = [""]
#predefined budgets
budget_shots = [
    [None],
    {50000, 1000, 333, 500, 250, 166, 125, 100, 83, 55, 41, 33, 27,
    500000, 10000, 3333, 5000, 2500, 1666, 1250, 1000, 833, 555, 416, 333, 277,
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
                    x // noise_levels
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
    """
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
    return storage


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

def clean_data(maxcut_df):
    """cleans up data and assigns noisy data label to appropiate data,
    also puts converged values in appropiate places. This is just the high
    budget ZNE and Noisy values that converged already to infinity at 10**5."""
    noisy = maxcut_df.query('abs_error > 0  & copies == 1 & nlsp==1 & res_type=="abs_error" & ( type=="VD")')
    noisy["type"] = "noisy"
    extra = maxcut_df.query('Qubits==4&budget==10**12&(type=="VD"|type=="ZNE")')
    extra_vd = maxcut_df.query('budget==10**12& type=="VD" & copies==1')
    extra_noise= maxcut_df.query('budget==10**12& type=="VD" & copies==1')
    extra_vd['budget']=10**10
    extra['budget']=10**10
    extra_noise['budget']=10**9
    return pd.concat([maxcut_df,extra,extra_vd,extra_noise])


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
        shot_budget = relative_shots(training, noise_levels, max_copies, budget)
        for deets in shot_budget.iterrows():
            copy_no = deets[1][0]
            technique = deets[1][1]
            allowance = deets[1][2]
            if budget == 0:
                df_temp_abs = df.query(f"budget=={budget} & copies == {copy_no}")
            else:
                df_temp_abs = df.query(
                    f'type=="{technique}" & shots == {allowance} &  copies == {copy_no}'
                ).assign(budget=budget)
            filtered_dfs.append(df_temp_abs)
    filtered_dfs = pd.concat(filtered_dfs)
    
    return filtered_dfs

def load_max_cut_data(qubit_number:int, training_set=100):
    from pathlib import Path
    if Path("./MaxCut_runs/full_data_MC{qubit_number}.pkl").is_file():
        Q=pd.read_pickle(f'./MaxCut_runs/full_data_MC{qubit_number}.pkl') 
    else:
        base_folder = f"./MaxCut_runs/Q{qubit_number}/"
        folders = [
            base_folder + ""
        ]  
        tags = [""]
        Q, _, _ = load_multiple_files_budget(
            [qubit_number],
            [0],
            29,
            nlsp_list=[1, 2, 3],
            N=10,
            Nts=training_set,
            max_copies=6,
            tags=tags,
            density_matrices=False,
            shots=budget_shots,
            folders=folders,
            train=True,
            budgets=[0, 1],
            train_use=100,
        )

        Q = Q.assign(description="3nlsp_full")
        Q.to_pickle(f'./MaxCut_runs/full_data_MC{qubit_number}.pkl')
    Q_filtered = filter_budget(Q,[0,10**5,10**6,10**7,10**8,10**9,10**10],100).drop(['shots'],axis=1).query("abs_error>0")
    
    ###
    Q_filtered.loc[Q_filtered["budget"] == 0, "budget"] = 10**12  
    ### Only done for plotting purposes, so that infinity is at the end of the plot rather than the start
    return Q_filtered

def figure_5(df,statistic_to_plot='mean'):
    """Plots figure five, which is the performance of the techniques over budget. This is extended in the sense that it returns
    all qubit results in a panel."""
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
        'abs_error > 0  &nlsp==1  & copies==4 & description == "3nlsp_full" & res_type=="abs_error" & ( type=="UNITED")'
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
        col_wrap = 2,
        estimator=statistic_to_plot,
        markers=True,
        ci=None,
    )
    fig.set(yscale="log", xscale="log")
    fig.set_ylabels(statistic_to_plot+' of absolute error')
    ax = fig.axes[0]
    xticks = ax.get_xticks().tolist()
    for i in range(len(xticks) - 1):
        xticks[i] = rf"$10^{{{i+3}}}$"
    xticks[-3] = r"$\infty$"
    xticks[-4] = r"$\cdots$"
    ax.set_xticklabels(xticks)
    return fig

def figure_4(maxcut_df,statistic_to_plot='mean',like_paper=True):
    """"Plot figure 5 from the paper, which is the absoulte error over qubits at different budgets."""
    if like_paper:
        df = maxcut_df.query('budget>0&budget<10**11&budget!=10**9&budget!=10**7')
    else: # Basically will show all computed budgets
        df = maxcut_df.query('budget>0&budget<10**11')
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
        'abs_error > 0  &nlsp==1  & copies==4& description == "3nlsp_full" & res_type=="abs_error" & ( type=="UNITED")'
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
        col='budget',
        col_wrap=2,
        y="abs_error",
        hue="type",
        estimator=statistic_to_plot,
        markers=True,
        ci=None,
    )
    fig.set_ylabels(statistic_to_plot+' of absolute error')
    fig.set(yscale="log")
    return fig

def load_all_data():
    results = {qubit:load_max_cut_data(qubit) for qubit in [4,6,8,10]}
    maxcut_df = pd.concat(results.values())
    maxcut_df = clean_data(maxcut_df)
    return maxcut_df

def load_unfiltered_data():
    dfs=[]
    for i in [4,6,8,10]:
        dfs.append(pd.read_pickle(f'./MaxCut_runs/full_data_MC{i}.pkl'))
    return pd.concat(dfs).reset_index()

def load_raw_maxcut_data():
    """Returns a dictionary of the coi data and training data with indices 'coi' and 'train'."""
    import os 
    
    def quick_load(fi, dir_path):
        res = []
        files = []
        for file in os.listdir(dir_path):
            if file.endswith('.pkl') and file.startswith(f"pandas_{fi}"):
                res.append(pd.read_pickle(dir_path+file).assign(file=file).assign(data=fi))
                files.append(file)
        res = pd.concat(res)
        res=res.reset_index()
        res = res.reset_index()
        res = res.drop(['result_type'],axis=1)
        res["exact"] = pd.to_numeric(res['exact'])
        res["expectation"] = pd.to_numeric(res['expectation'])
        res["exact_abs"] = np.abs(res['exact'])
        return res
    
    train = []
    coi = []
    for i in [4,6,8,10]:
        path = f'./MaxCut_runs/Q{i}/'
        train.append(quick_load('train',path))
        coi.append(quick_load('COI',path))
    train = pd.concat(train).reset_index().drop('level_0',axis=1)
    coi = pd.concat(coi).reset_index().drop('level_0',axis=1)
    return {'coi':coi,'train':train}

def figure_7(df, budget=10**10, statistic_to_plot='mean'):
    """Plots figure 7 for any budget"""
    df = df.query(f'budget=={budget}')
    few_copy_methods = df.query(
        'abs_error >  0  & nlsp==1  & description == "3nlsp_full" & res_type=="abs_error" & ( type=="VD")'
    )
    many_copy_methods = df.query(
        'abs_error > 0  &nlsp==1  & description == "3nlsp_full" & res_type=="abs_error" & ( type=="UNITED")'
    )
    plot_df = pd.concat(
        [few_copy_methods, many_copy_methods],
        axis=0,
        ignore_index=True,
    )
    fig = sns.relplot(
        data=plot_df.reset_index(),
        kind="line",
        x="copies",
        y="abs_error",
        hue="Qubits",
        style='type',
        estimator=statistic_to_plot,
        ci=None,
    )
    fig.set(yscale='log',title=f'budget {budget}')
    fig.set_ylabels(statistic_to_plot+' of absolute error')
    return fig
    
def figure_8(dataset):
    """Plots the disctibution of training data and circuit of interest for desired qubit"""
    sns.displot(kind='hist',data=dataset,x='exact',col='qubits',col_wrap=2,bins=25,facet_kws=dict(sharey=True)).set(yscale='log')

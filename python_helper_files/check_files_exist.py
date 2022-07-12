# %%
import os

test_budgets = [
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
]
qubits = [4, 4, 6, 6, 8, 8, 10]
depths = [4, 4 * 16, 6, 6 * 16, 8, 8 * 16, 10]
nlsps = [1, 2, 3]
N = 10
Nts = 100
tag = "_LinGrow"
seed = range(30)


def missing_file_checker(
    Q=qubits,
    p=depths,
    num_seeds=30,
    nlsp_list=[1, 2, 3],
    N=10,
    Nts=100,
    max_copies=6,
    tags=["_LinGrow"],
    shots=test_budgets,
    folder="./results/",
    train=True,
    budgets=[0, 10 ** 5, 10 ** 6, 10 ** 7, 10 ** 8],
):
    """Generates a list of missing files per budget.

    Args:
        Q (int, optional): qubits being looked at. Defaults to qubits.
        p (int, optional): deoth of files. Defaults to depths.
        num_seeds (int, optional): number of seeds checked. Defaults to 30.
        nlsp_list (list, optional): list of noise levels. Defaults to [1,2,3].
        N (int, optional): number of non cliffords. Defaults to 10.
        Nts (int, optional): number of training samples. Defaults to 100.
        max_copies (int, optional): copies for VD. Defaults to 6.
        tags (list, optional): additional tags. Defaults to ['_LinGrow'].
        shots (int, optional): shots in a given allowance. Defaults to test_budgets.
        folder (str, optional): folder. Defaults to './results/'.
        train (bool, optional): include training data?. Defaults to True.
        budgets (list, optional): name of budgets, can be arbitrary. Defaults to [0, 10**5,10**6,10**7,10**8].

    Returns:
        list: list of list with strings of missing files.
    """
    missing_file_list = [[], [], [], [], []]
    file = folder + "ALL.txt"
    with open(file, "w") as f:
        for budget_i, budget in enumerate(budgets):
            for shot in shots[budget_i]:
                for qubit_no, depth in zip(Q, p):
                    for nlsp in nlsp_list:
                        for seed_i in range(num_seeds):
                            for tag in tags:
                                if shot != None:
                                    train = (
                                        folder
                                        + "train_data_Q{}p{}N{}Nts{}MC{}nlsp{}seed{}shots{}{}.mat".format(
                                            qubit_no,
                                            depth,
                                            N,
                                            Nts,
                                            max_copies,
                                            nlsp,
                                            seed_i,
                                            shot,
                                            tag,
                                        )
                                    )
                                    follow = (
                                        folder
                                        + "coi_data_Q{}p{}MC{}nlsp{}seed{}shots{}{}.mat".format(
                                            qubit_no,
                                            depth,
                                            max_copies,
                                            nlsp,
                                            seed_i,
                                            shot,
                                            tag,
                                        )
                                    )
                                    namtrain = "train_data_Q{}p{}N{}Nts{}MC{}nlsp{}seed{}shots{}{}.mat".format(
                                        qubit_no,
                                        depth,
                                        N,
                                        Nts,
                                        max_copies,
                                        nlsp,
                                        seed_i,
                                        shot,
                                        tag,
                                    )
                                    namfollow = "coi_data_Q{}p{}MC{}nlsp{}seed{}shots{}{}.mat".format(
                                        qubit_no,
                                        depth,
                                        max_copies,
                                        nlsp,
                                        seed_i,
                                        shot,
                                        tag,
                                    )

                                else:
                                    train = (
                                        folder
                                        + "train_data_Q{}p{}N{}Nts{}MC{}nlsp{}seed{}{}.mat".format(
                                            qubit_no,
                                            depth,
                                            N,
                                            Nts,
                                            max_copies,
                                            nlsp,
                                            seed_i,
                                            tag,
                                        )
                                    )
                                    follow = (
                                        folder
                                        + "coi_data_Q{}p{}MC{}nlsp{}seed{}{}.mat".format(
                                            qubit_no,
                                            depth,
                                            max_copies,
                                            nlsp,
                                            seed_i,
                                            tag,
                                        )
                                    )
                                    namtrain = "train_data_Q{}p{}N{}Nts{}MC{}nlsp{}seed{}shots{}{}.mat".format(
                                        qubit_no,
                                        depth,
                                        N,
                                        Nts,
                                        max_copies,
                                        nlsp,
                                        seed_i,
                                        shot,
                                        tag,
                                    )
                                    namfollow = "coi_data_Q{}p{}MC{}nlsp{}seed{}shots{}{}.mat".format(
                                        qubit_no,
                                        depth,
                                        max_copies,
                                        nlsp,
                                        seed_i,
                                        shot,
                                        tag,
                                    )

                                if not os.path.isfile(train):
                                    missing_file_list[budget_i].append(train)
                                    f.write(namtrain + "\n")
                                if not os.path.isfile(follow):
                                    missing_file_list[budget_i].append(follow)
                                    f.write(namfollow + "\n")
    return missing_file_list


# %%

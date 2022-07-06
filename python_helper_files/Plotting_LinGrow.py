# %%
import enum
from math import gamma
import scipy.io as sio
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt6
from scipy.io import savemat
import random
from seaborn import palettes 
from sklearn import linear_model
import pandas as pd
import os
import qutip as q
import matplotlib.pyplot as plt
plt.style.use(['science','no-latex'])
from post_processing_multi import * 

import seaborn as sns
sns.set_context("paper")
def relative_shots(training_sets,noise_levels,max_copies,shots_budget):
    if shots_budget >= 10**9:
        CDR_shots =  shots_budget//training_sets
        vnCDR_shots = CDR_shots//noise_levels
        CGVD_shots =  shots_budget//training_sets
        storage = pd.DataFrame({
            "VD":[shots_budget]+[shots_budget//2]*(max_copies-1),
            "ZNE":[shots_budget//2]+[shots_budget//4]*(max_copies-1),
            "CDR":[CDR_shots]+[CDR_shots//2]*(max_copies-1),
            "vnCDR":[vnCDR_shots]*(max_copies),
            "CGVD":[base_cost//(copy+1)//2 for copy,base_cost in enumerate([CGVD_shots]*(max_copies))],
            "vnCGVD":[x//noise_levels//2 for x in [base_cost//(copy+1) for copy,base_cost in enumerate([CGVD_shots]*(max_copies))]],
            "copies":list(range(1,max_copies+1))
        })
    elif shots_budget != 0:
        CDR_shots =  shots_budget//training_sets
        vnCDR_shots = CDR_shots//noise_levels
        CGVD_shots =  shots_budget//training_sets
        storage = pd.DataFrame({
            "VD":[shots_budget//2]*(max_copies),
            "ZNE":[shots_budget//2]+[shots_budget//4]*(max_copies-1),
            "CDR":[CDR_shots]+[CDR_shots//2]*(max_copies-1),
            "vnCDR":[vnCDR_shots]*(max_copies),
            "CGVD":[base_cost//(copy+1)//2 for copy,base_cost in enumerate([CGVD_shots]*(max_copies))],
            "vnCGVD":[x//noise_levels//2 for x in [base_cost//(copy+1) for copy,base_cost in enumerate([CGVD_shots]*(max_copies))]],
            "copies":list(range(1,max_copies+1))
        })

    else:        
        storage = pd.DataFrame({
            "VD":None,
            "ZNE":None,
            "CDR":[None]*6,
            "vnCDR":[None]*6,
            "CGVD":[None]*6,
            "vnCGVD":[None]*6,
            "copies":list(range(1,max_copies+1))})

        
    return pd.melt(storage,id_vars=('copies'))
def relative_shots_table(training_sets,noise_levels,max_copies,shots_budget,swap_levels):
    if shots_budget != 0:
        CDR_shots = shots_budget//training_sets
        vnCDR_shots = CDR_shots//noise_levels
        CGVD_shots = shots_budget//training_sets
        storage = pd.DataFrame({
            "VD":[shots_budget]+[shots_budget//2]*(max_copies-1),
            "ZNE":[shots_budget//2]+[shots_budget//4]*(max_copies-1),
            "CDR":[CDR_shots]+[CDR_shots//2]*(max_copies-1),
            "vnCDR":[vnCDR_shots]+[vnCDR_shots//2]*(max_copies-1),
            "CGVD":[CGVD_shots]+[base_cost//(copy+2)//2 for copy,base_cost in enumerate([CGVD_shots]*(max_copies-1))],
            "vnCGVD":[vnCDR_shots]+[x//noise_levels//2 for x in [base_cost//(copy+2) for copy,base_cost in enumerate([CGVD_shots]*(max_copies-1))]],
            "vnCGVD+":[vnCDR_shots//swap_levels]+[x//noise_levels//2//swap_levels for x in [base_cost//(copy+2) for copy,base_cost in enumerate([CGVD_shots]*(max_copies-1))]],
            "copies":list(range(1,max_copies+1))
        })
    else:
        storage = pd.DataFrame({
            "VD":None,
            "ZNE":None,
            "CDR":None,
            "vnCDR":None,
            "CGVD":None,
            "vnCGVD":None,
            "copies":list(range(1,max_copies+1))})

    return storage
def cost_fucntion_eval_table(training_sets,noise_levels,max_copies,shots_budget):
    CDR_shots = shots_budget//training_sets
    vnCDR_shots = CDR_shots//noise_levels
    CGVD_shots = shots_budget//training_sets
    storage = pd.DataFrame({
        "VD":[1]+[2]*(max_copies-1),
        "ZNE":[2]+[4]*(max_copies-1),
        "CDR":[100]+[200]*(max_copies-1),
        "vnCDR":[300]+[600]*(max_copies-1),
        "CGVD":[100]+[base_cost*(copy+2)*2 for copy,base_cost in enumerate([100]*(max_copies-1))],
        "vnCGVD":[300]+[x*noise_levels*2 for x in [base_cost*(copy+2) for copy,base_cost in enumerate([100]*(max_copies-1))]],
        "copies":list(range(1,max_copies+1))
    })
    return storage

def single_qubit_gate_length(layers,qubits):
    return layers*6*(qubits-1)

def two_qubit_gate_length(layers,qubits):
    e2=layers*(qubits-1)
    
def circuit_length(layers,qubits):
    e1=layers*6*(qubits-1)
    e2=layers*(qubits-1)
    return e1+e2

def total_error_rates_normal(layers,qubits,factor=1):
    e1=layers*6*(qubits-1)*0.0011*factor
    e2=layers*(qubits-1)*0.0021*factor
    return e1+e2

def load_multiple_files(Q,p,num_seeds,nlsp_list = [1,2,3],N=10,Nts=100,max_copies=6,tags=[''],
                               density_matrices = False,shots=[None],folder='',train=True):
    dfs_error_rescaling = []
    dfs_absolute_error  = [] 
    dfs_standard_values = []
    for qubit_no,depth in zip(Q,p):
        print(f"Qubit: {Q}, depth: {p}")
        data_processing =  VD_CDR(
                            qubit_no,depth,num_seeds,N,Nts,max_copies,nlsp_list,shots=shots,
                            folder=folder,autoprocessing=True,plots=False,
                            density_matrices = density_matrices,extra_tags=tags,train=train)
                        
        dfs_absolute_error.append(
            data_processing.abs_error_df.assign(Qubits= qubit_no, depth= depth)
            )
        dfs_error_rescaling.append(
            data_processing.calculated_vals_df.assign(Qubits= qubit_no, depth= depth)
            )
        dfs_standard_values.append(
            data_processing.standard_df.assign(Qubits= qubit_no, depth= depth)
            )
    dfs_absolute_error= pd.concat(dfs_absolute_error)
    dfs_error_rescaling= pd.concat(dfs_error_rescaling)
    dfs_standard_values= pd.concat(dfs_standard_values)

    return dfs_absolute_error, dfs_error_rescaling, dfs_standard_values

def load_multiple_files_budget(Q,p,num_seeds,nlsp_list = [1,2,3],N=10,Nts=100,max_copies=6,tags=[''],
                               density_matrices = False,shots=[[None]],folders=[''],train=True,budgets=[1],train_use=100):
    dfs_error_rescaling = []
    dfs_absolute_error  = [] 
    dfs_standard_values = []
    for tag,folder in zip(tags,folders):
        print(tag)
        for budget_i, budget in enumerate(budgets):
            for qubit_no,depth in zip(Q,p):
                    print(f"Qubit: {qubit_no}, depth: {depth}, budget:{budget}")
                    # try:
                    data_processing =  VD_CDR(
                                        qubit_no,depth,num_seeds,N,Nts,max_copies,nlsp_list,shots=shots[budget_i],
                                        folder=folder,autoprocessing=True,
                                        density_matrices = density_matrices,extra_tags=[tag],train=train,budget=budget,train_use=train_use)                                    
                    dfs_absolute_error.append(
                        data_processing.abs_error_df.assign(Qubits= qubit_no, depth= depth,budget=budget)
                        )
                    dfs_error_rescaling.append(
                        data_processing.calculated_vals_df.assign(Qubits= qubit_no, depth= depth,budget=budget)
                        )
                    dfs_standard_values.append(
                        data_processing.standard_df.assign(Qubits= qubit_no, depth= depth,budget=budget)
                        )
                    # except: 
                    #     print("FAILED")
                    #     pass
    dfs_absolute_error = pd.concat(dfs_absolute_error)
    dfs_error_rescaling= pd.concat(dfs_error_rescaling)
    dfs_standard_values= pd.concat(dfs_standard_values)
    # dfs_absolute_error.to_pickle( folder+f"abs_error.pkl")
    # dfs_error_rescaling.to_pickle(folder+f"error_rescaling.pkl")
    # dfs_standard_values.to_pickle(folder+f"raw_results.pkl")

    return dfs_absolute_error, dfs_error_rescaling, dfs_standard_values
def filter_budget(df,budgets):
    filtered_dfs = []
    for budget_i, budget in enumerate(budgets):
        shot_budget = relative_shots(100,3,6,budget)
        for deets in shot_budget.iterrows():
            copy_no   = deets[1][0]
            technique = deets[1][1]
            allowance = deets[1][2]
            if budget ==0:
                df_temp_abs = df.query(f'budget=={budget} & type.str.startswith("{technique}") & copies == {copy_no}')
            else:
                df_temp_abs = df.query(f'shots == {allowance} & type.str.startswith("{technique}") & copies == {copy_no}').assign(budget=budget)
            filtered_dfs.append(df_temp_abs)

    filtered_dfs = pd.concat(filtered_dfs)
    return filtered_dfs

# %% 
# budge_abs = pd.read_pickle('./new_data/BIG_CONGLOMERATE_SHOTS/BIG_CONGLOMERATE/'+"abs_error.pkl")

base_folder = 'D:/Databases/VD_CDR/BIG_CONGLOMERATE/'
folders = [base_folder+'/LinGrow_complete/realistic_noise_model_shots/']#,base_folder+'LinGrow_complete/realistic_noise_model_shots/']#,base_folder+'RNM_e0p1/',base_folder+'RNM_e1/']
tags = ['_LinGrow']
budgies = [ [None],
            [50000,  1000,   333,   500,   250,   166,   125,   100,    83, 55,    41,    33,    27],
             [500000,  10000,   3333,   5000,   2500,   1666,   1250,   1000, 833,    555,    416,    333,    277],
             [5000000,  100000,   33333,   50000,   25000,   16666,   12500, 10000,    8333,    5555,    4166,    3333,    2777] ,
             [50000000,  1000000,   333333,   500000,   250000,   166666, 125000,   100000,    83333,    55555,    41666,    33333, 27777],
             [500000000,  10000000,   3333333,   5000000,   2500000,   1666666, 1250000,   1000000,    833333,    555555,    416666,    333333, 277777],
             [10000000000, 5000000000,  100000000,   33333333,   50000000,   25000000,   16666666, 12500000,   10000000,    8333333,    5555555,    4166666,    3333333, 2777777],
             [50000000000,  1000000000,   333333333,   500000000,   250000000,   166666666, 125000000,   100000000,    83333333,    55555555,    41666666,    33333333, 27777777],
             ]
big_budge = [[1000000000,   333333333,   250000000,   166666666, 125000000,   83333333,    55555555,    41666666,    27777777]]
budge_abs,budge_res,budge_val = load_multiple_files_budget([4,4,6,6,8,8,10],[4,4*16,6,6*16,8,8*16,10],30,nlsp_list = [1,2,3],N=10,Nts=100,max_copies=6,tags=tags,
                               density_matrices = False,shots=budgies,folders = folders,train=True,budgets=[0,5,6,7,8,9,10,11],train_use=100) 
NLSP3FULL_abs = budge_abs
NLSP3HALF_abs  ,_,_ = load_multiple_files_budget([4,4,6,6,8,8,10],[4,4*16,6,6*16,8,8*16,10],30,nlsp_list = [1,2,3],N=10,Nts=100,max_copies=6,tags=tags,
                               density_matrices = False,shots=budgies,folders = folders,train=True,budgets=[0,5,6,7,8,9,10,11],train_use=50) 
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
NLSP3FULL_abs  = NLSP3FULL_abs.assign(description =     '3nlsp_full')
# NLSP3QUART_abs = NLSP3QUART_abs.assign(description = '3nlsp_quar')
NLSP3HALF_abs  = NLSP3HALF_abs.assign(description =  '3nlsp_half')

# NLSP2HALF_abs ,
# NLSP2FULL_abs ,
# NLSP2QUART_abs,
# NLSP3QUART_abs,
# ,
everything= pd.concat([NLSP3FULL_abs,NLSP3HALF_abs])

# %%
everything_full = everything
everything_full.to_pickle(base_folder+"EVERYTHING_FULL_more_Data_complete_intercept.pkl")
# %%
# budge_abs = pd.read_pickle('./BIG_CONGLOMERATE/'+"LINGROW2nlsp_abs_0filteredfull_error.pkl")
# everything_full = pd.read_pickle(base_folder+"EVERYTHING_FULL_more_Data_filt.pkl")

t = filter_budget(everything_full,[0,10**5,10**6,10**7,10**8,10**9,10**10,10**11])
t['volume'] = t.depth*t.Qubits*t.Qubits
t["g"] = t.depth//t.Qubits
# %%

t= pd.read_pickle(base_folder+"PAPER_DATA_v3.pkl")

# t = pd.read_pickle('./BIG_CONGLOMERATE/'+"EVERYTHING_no_swap_nosie_filtered.pkl")
# %%
scaled_SWAP_Results = []
for i in [0]:
    for QUBIT,g in zip([4],[5]):
        temp = budge_abs.query(f'abs_error != 0 & nlsp == 1 & budget == {i} & Qubits == {QUBIT} & depth == {g}').copy()
        temp['scaled_error'] = temp['abs_error']/budge_abs.query(f'abs_error != 0 & nlsp == 1 & budget == {i} & copies == 1 & type == "VD_LinGrow" & Qubits == {QUBIT} & depth == {g}')['abs_error'].mean()
        scaled_SWAP_Results.append(temp)
scaled_swap_panda = pd.concat(scaled_SWAP_Results)

# %%
scaled_Results = []
for i in [0,10**5,10**6,10**7,10**8,10**9,10**10,10**11]:
    for description in ['3nlsp_full','3nlsp_half']:
        for QUBIT,g in zip([4,4,6,6,8,8,10],[1,16,1,16,1,16,1]):
            temp = t.query(f'abs_error != 0 & description == "{description}" & nlsp == 1 & budget == {i} & Qubits == {QUBIT} & g == {g}').copy()
            temp['scaled_error'] = temp['abs_error']/t.query(f'abs_error != 0  & description == "{description}" & nlsp == 1 & budget == {i} & copies == 1 & type == "VD_LinGrow" & Qubits == {QUBIT} & g == {g}')['abs_error'].mean()
            scaled_Results.append(temp)
scaled_panda = pd.concat(scaled_Results)

# scaled_panda.to_pickle(base_folder+"EVERYTHING_FULL_more_files_filtered.pkl")
# scaled_Results = []
# for i in [0,10**5,10**6,10**7,10**8,10**9,10**10,10**11]:
#     for description in ['3nlsp_full','3nlsp_half']:
#         for QUBIT,g in zip([4,4,6,6,8,8,10],[1,16,1,16,1,16,1]):
#             temp = t.query(f'abs_error != 0 &  nlsp == 1 & budget == {i} & Qubits == {QUBIT} & g == {g}').copy()
#             temp['scaled_error'] = temp['abs_error']/t.query(f'abs_error != 0  &  nlsp == 1 & budget == {i} & copies == 1 & type == "VD_LinGrow" & Qubits == {QUBIT} & g == {g}')['abs_error'].mean()
#             scaled_Results.append(temp)
# scaled_panda = pd.concat(scaled_Results)
# scaled_panda.to_pickle(base_folder+"EVERYTHING_FULL_more_files_filtered.pkl")
# %%
scaled_panda = pd.read_pickle(base_folder+"EVERYTHING_FULL_more_files_filtered.pkl")
# %%
base_folder = 'D:/Databases/VD_CDR/'

t['type'] = t['type'].str.replace(r'_LinGrow', '')

ticker =t.type.unique().tolist()+['Noisy']
colors = sns.color_palette('tab10', n_colors=len(ticker))  # get a number of colors
cmap = dict(zip(ticker, colors))  # zip values to colors
marker_map = {'VD':'s','ZNE':'*','CDR':'h','vnCDR':'^','vnCGVD':'p','Noisy':'o','CGVD':'v'}
cmap = {'VD':'green','ZNE':'red','CDR':'pink','vnCDR':'purple','vnCGVD':'blue','Noisy':'black','CGVD':'orange'}
line_map = {'VD':(4, 1.5),'ZNE':(4, 1.5),'CDR':(4, 1.5),'vnCDR':(4, 1.5),'vnCGVD':(4, 1.5),'Noisy':'','CGVD':(4, 1.5)}

# %%
tags = ['_LinGrow']
title_dict = {'_LinGrow':'Linear Growth','_RNM_e1':'Scaled Noise (e_r=1)','_RNM_e0p1':'Scaled Noise (e_r=0.1)'}
for tag in tags:
    for QUBIT,g in zip([4,4,6,6,8,8,10],[1,16,1,16,1,16,1]):
        fig, axes = plt.subplots(2,3,figsize=(15,10),sharey='row',sharex=True)
        axes=axes.flatten()
        dfs = []
        for i in [10**7,10**8,10**9,10**10]:
            df1 = scaled_panda.query(
            f'budget=={i}    & Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0 & description.str.startswith("3") & description.str.endswith("full")')
            df2 = scaled_panda.query(
            f'budget=={i*10}    & Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0 & ~(description.str.startswith("3") & description.str.endswith("full"))')
            dfs.append(pd.concat([df1,df2]))
        sns.lineplot(data=scaled_panda.query(
            f'budget==0 &   Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0      &  type.str.endswith("{tag}")& type.str.startswith("vn")'),
            ax=axes[0],x='copies',y='scaled_error', style='type',hue='description' ,estimator='mean',markers=marker_map,err_style=None,palette=cmap)
        sns.lineplot(data=scaled_panda.query(
            f'budget==0 &   Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0      &  type.str.endswith("{tag}")& type.str.startswith("vn")'),
            ax=axes[3],x='copies',y='scaled_error', style='type',hue='description' ,estimator='mean',markers=marker_map,err_style=None,palette=cmap)

        sns.lineplot(data=dfs[0].query(
            f'   Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0      &  type.str.endswith("{tag}")& type.str.startswith("vn")'),
            ax=axes[0+1],x='copies',y='scaled_error', style='type',hue='description' ,estimator='mean',markers=marker_map,err_style=None,palette=cmap)
        sns.lineplot(data=dfs[1].query(
            f'    Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0 &  type.str.endswith("{tag}") & type.str.startswith("vn")'),
            ax=axes[1+1],x='copies',y='scaled_error', style='type',hue='description' ,estimator='mean',markers=marker_map,err_style=None,palette=cmap)
        sns.lineplot(data=dfs[2].query(
            f' Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0 &  type.str.endswith("{tag}") & type.str.startswith("vn")'),
            ax=axes[2+2],x='copies',y='scaled_error', style='type',hue='description' ,estimator='mean',markers=marker_map,err_style=None,palette=cmap)
        sns.lineplot(data=dfs[3].query(
            f' Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0 &  type.str.endswith("{tag}") & type.str.startswith("vn")'),
            ax=axes[3+2],x='copies',y='scaled_error', style='type',hue='description' ,estimator='mean',markers=marker_map,err_style=None,palette=cmap)

        # sns.lineplot(data=scaled_panda.query(
            # f'budget==1000000   & Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0 &  type.str.endswith("{tag}")'),
            # ax=axes[2],x='copies',y='scaled_error', hue='type',style='type' ,estimator='mean',err_style=None)
        # sns.lineplot(data=scaled_panda.query(
            # f'budget==10000000  & Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0 &  type.str.endswith("{tag}")'),
            # ax=axes[3],x='copies',y='scaled_error', hue='type',style='type' ,estimator='mean',err_style=None)

        for i in axes:
            i.set_yscale('log')
            i.set_ylim([10**(-2),10**1])
            i.legend([],[] ,frameon=False)

        axes[1].set_title('10^7   Budget')
        axes[2].set_title('10^8   Budget')
        axes[4].set_title('10^9   Budget')
        axes[5].set_title('10^10   Budget')

        handles, labels = axes[0].get_legend_handles_labels()

        fig.legend(handles, labels, title='Mitigation Strategies:',loc='lower center', bbox_to_anchor=(0.5, -0.025),
                ncol=7, fancybox=True, shadow=True)

        # axes[2].set_title('1M  Budget')
        # axes[3].set_title('10M Budget')
        fig.suptitle(f"different noise level and training data at {title_dict[tag]} Q={QUBIT} g={g}")
        plt.savefig(base_folder+f'/PLOTS/MT1_nlsp{tag}_budgetsQ{QUBIT}g{g}.png')
        plt.close(fig)
 # %%
o = t.query(
    '(budget==10**5|budget==10**10) & res_type=="abs_error" & nlsp == 1 & abs_error != 0   & description == "3nlsp_full"  & type.str.startswith("VD")')

g = sns.FacetGrid(o, row="budget", sharey='row')
g.map_dataframe(sns.lineplot, x="copies", y="abs_error", hue="volume",style="volume",err_style=None,palette='tab10').set(yscale ='log')
g.add_legend(bbox_to_anchor=(0.85, 0.3))
g.set_axis_labels("Copies", r"Mean $\langle \sigma^1_z \rangle$ absolute error")

 
 # %%
o = t.query(
    'res_type=="abs_error" & nlsp == 1 & abs_error != 0  & budget != 0 description == "3nlsp_full" & copies == 1 & type.str.startswith("VD")').copy()
o['type'] = 'Noisy'
a = t.query(
    'res_type=="abs_error" & nlsp == 1 & abs_error != 0  & budget != 0 & description == "3nlsp_full" & copies == 1 & type.str.startswith("ZNE")').copy().reset_index(drop=True)
b = t.query(
    '
     res_type=="abs_error" &  nlsp == 1 & abs_error != 0  & budget != 0 & description == "3nlsp_full" & copies == 2 & type.str.startswith("VD")').copy()
c = t.query(
    'res_type=="abs_error" &  nlsp == 1 & abs_error != 0  & budget != 0 & description == "3nlsp_half" & copies == 3 & type.str.startswith("vnCG")').copy()
d = t.query(
    'res_type=="abs_error" & nlsp == 1 & abs_error != 0  & budget != 0 & description == "3nlsp_half" & copies == 1 & type.str.startswith("vnCD")').copy()
d['shots']=2*d['shots']
c['shots']=2*c['shots']

e = t.query(
    'g==16 & Qubits==4 & res_type=="abs_error" &  nlsp == 1 & abs_error != 0  & budget != 0 & description == "3nlsp_full" & copies == 4 & type.str.startswith("VD")').copy()
f = t.query(
    'g==16 & Qubits!=4 & res_type=="abs_error" &  nlsp == 1 & abs_error != 0  & budget != 0 & description == "3nlsp_full" & copies == 2 & type.str.startswith("VD")').copy()


idnex_mis = a.query('budget==10**9 & Qubits == 6 & g == 1').index
idnex_rep = a.query('budget==10**8 & Qubits == 6 & g == 1').index
a.loc[idnex_mis,'abs_error'] = list(a.iloc[idnex_rep]['abs_error'])
idnex_mis = a.query('budget==10**9 & Qubits == 10 & g == 1').index
idnex_rep = a.query('budget==10**8 & Qubits == 10 & g == 1').index
a.loc[idnex_mis,'abs_error'] = list(a.iloc[idnex_rep]['abs_error'])

c['budget'] = c['budget']//10
d['budget'] = d['budget']//10
c['shots'] = c['shots']//10
d['shots'] = d['shots']//10

selection = pd.concat([a,b,c,d,o,e,f])
selection['type'] = selection['type'].str.replace(r'_LinGrow', '')
 
# %%
o = t.query(
    ' nlsp == 1 & abs_error != 0   & description == "3nlsp_full"  & type.str.startswith("VD")').copy()
o['type'] = 'Noisy'
a = t.query(
    ' nlsp == 1 & abs_error != 0   & description == "3nlsp_full"  & type.str.startswith("ZNE")').copy().reset_index(drop=True)
b = t.query(
    ' abs_error != 0   & description == "3nlsp_full"  & type.str.startswith("VD")').copy()
c = t.query(
    '  nlsp == 1 & abs_error != 0   & description == "3nlsp_half"  & type.str.startswith("vnCG")').copy()
d = t.query(
    ' nlsp == 1 & abs_error != 0   & description == "3nlsp_half"  & type.str.startswith("vnCD")').copy()
i = t.query(
    ' abs_error != 0   & description == "3nlsp_half"  & type.str.startswith("CG")').copy()
j = t.query(
    ' abs_error != 0   & description == "3nlsp_half"  & type.str.startswith("CD")').copy()
d['shots']=2*d['shots']
c['shots']=2*c['shots']
i['shots']=2*i['shots']
j['shots']=2*j['shots']
e = t.query(
    '  nlsp == 1 & abs_error != 0   & description == "3nlsp_full"  & type.str.startswith("vnCG")').copy()
f = t.query(
    ' nlsp == 1 & abs_error != 0   & description == "3nlsp_full"  & type.str.startswith("vnCD")').copy()
g = t.query(
    '  abs_error != 0   & description == "3nlsp_full"  & type.str.startswith("CG")').copy()
h = t.query(
    ' abs_error != 0   & description == "3nlsp_full"  & type.str.startswith("CD")').copy()

idnex_mis = a.query('copies == 1  & budget==10**9 & Qubits == 6 & g == 1').index
idnex_rep = a.query('copies == 1  & budget==10**8 & Qubits == 6 & g == 1').index
a.loc[idnex_mis,'abs_error'] = list(a.iloc[idnex_rep]['abs_error'])
idnex_mis = a.query('copies == 1  & budget==10**9 & Qubits == 10 & g == 1').index
idnex_rep = a.query('copies == 1  & budget==10**8 & Qubits == 10 & g == 1').index
a.loc[idnex_mis,'abs_error'] = list(a.iloc[idnex_rep]['abs_error'])

c['budget'] = c['budget']//10
d['budget'] = d['budget']//10
c['shots'] = c['shots']//10
d['shots'] = d['shots']//10
j['budget'] = j['budget']//10
i['budget'] = i['budget']//10
i['shots'] = i['shots']//10
j['shots'] = j['shots']//10

selection = pd.concat([a,b,c,d,o,e,f,g,h,i,j])
selection['type'] = selection['type'].str.replace(r'_LinGrow', '')
selection['type'] = selection['type'].str.replace('vnCGVD', 'UNITED')

# %%

g = sns.FacetGrid(selection, col="Qubits", sharey='row',  row="g")
g.map_dataframe(sns.lineplot, x="budget", y="abs_error", hue="type",style="type",dashes=line_map,markers=marker_map,err_style=None,palette=cmap).set(yscale ='log', xscale ='log',xlim=[10**5,10**10+10**9],xticks=[10**5,10**6,10**7,10**8,10**9,10**10])
g.add_legend(bbox_to_anchor=(0.85, 0.3))
g.set_axis_labels("Budget", r"Mean $\langle \sigma^1_z \rangle$ absolute error")
g.fig.delaxes(g.axes[-1, -1])
g.axes[0][0].set_ylabel('')
g.axes[1][0].yaxis.set_label_coords(-0.125,1.08)
for ax in g.axes[0]:
    ax.set_ylim([10**-4,10**-1])
for ax in g.axes[1]:
    ax.set_ylim([10**-2,10**0])
plt.tight_layout()
# plt.savefig(base_folder+'abs_error_per_budget.pdf')
# %%



o_ = t.query(
    'res_type=="abs_error" & nlsp == 1 & abs_error != 0  & description == "3nlsp_full" & copies == 1 & type.str.startswith("VD")').copy()
o_['type'] = 'Noisy'
a_ = t.query(
    'res_type=="abs_error" & nlsp == 1 & abs_error != 0  & description == "3nlsp_full" & copies == 1 & type.str.startswith("ZNE")').copy()
b_ = t.query(
    'res_type=="abs_error" & nlsp == 1 & abs_error != 0  & description == "3nlsp_full" & copies == 2 & type.str.startswith("VD")').copy()
c_ = t.query(
    'res_type=="abs_error" & nlsp == 1 & abs_error != 0  & description == "3nlsp_half" & copies == 3 & type.str.startswith("vnCG")').copy()
d_ = t.query(
    'res_type=="abs_error" & nlsp == 1 & abs_error != 0  & description == "3nlsp_half" & copies == 1 & type.str.startswith("vnCD")').copy()

c_['budget'] = c_['budget']//10
d_['budget'] = d_['budget']//10
c_['shots'] = c_['shots']*2
d_['shots'] = d_['shots']*2

o_['budget'] = np.log10(o_['budget'])
c_['budget'] = np.log10(c_['budget'])
d_['budget'] = np.log10(d_['budget'])
a_['budget'] = np.log10(a_['budget'])
b_['budget'] = np.log10(b_['budget'])



o = t.query(
    'res_type=="abs_error" & nlsp == 1 & abs_error != 0 & description == "3nlsp_full" & copies == 1 & type.str.startswith("VD")').copy()
o['type'] = 'Noisy'
a = t.query(
    'res_type=="abs_error" & nlsp == 1 & abs_error != 0 & description == "3nlsp_full" & copies == 1 & type.str.startswith("ZNE")').copy()
b = t.query(
    'res_type=="abs_error" & nlsp == 1 & abs_error != 0 & description == "3nlsp_full" & copies == 2 & type.str.startswith("VD")').copy()
c = t.query(
    'res_type=="abs_error" & nlsp == 1 & abs_error != 0  & description == "3nlsp_half" & copies == 3 & type.str.startswith("vnCG")').copy()
d = t.query(
    'res_type=="abs_error" & nlsp == 1 & abs_error != 0  & description == "3nlsp_half" & copies == 1 & type.str.startswith("vnCD")').copy()

c['budget'] = c['budget']//10
d['budget'] = d['budget']//10
c['shots'] = c['shots']*2
d['shots'] = d['shots']*2

o['budget'] = np.log10(o['budget'])
c['budget'] = np.log10(c['budget'])
d['budget'] = np.log10(d['budget'])
a['budget'] = np.log10(a['budget'])
b['budget'] = np.log10(b['budget'])

selection = pd.concat([a,b,c,d,o,a_,b_,c_,d_,o_])
# %% 
selection = pd.read_csv("dataframe_MEAN_VD_COPIES_PLOTTING.txt")
g = sns.FacetGrid(selection.query("result_type=='abs_error'"), col="budget")
g.map_dataframe(sns.lineplot, x="copies", y="value", hue="qubit",style='g').set(yscale ='log')
g.add_legend(bbox_to_anchor=(0.85, 0.3))
g.set_axis_labels("Budget", r"Mean $\langle \sigma^1_z \rangle$ absolute error")
# plt.savefig(base_folder+'abs_error_over_budget.pdf')
 # %%
import csv
with open('dataframe_MEAN_ALT_COPIES_PLOTTING.txt', 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter=',')
    line_count = 0
    csv_writer.writerow(['type', 'budget', 'qubit', 'g','value','shots','copies','description','result_type','id_string','volume'])
    for thing in selection.type.unique():
        print(thing)
        for copies in [1]:
            for budget in [10**5,10**10]:
                for q, g in zip([4,6,8,10,4],[1,1,1,1,16]):
                    for result_type in ['abs_error']:
                        x = selection.query(f" res_type == '{result_type}' &budget == {budget} & type == '{thing}' & Qubits == {q} & g == {g}")['abs_error'].mean()
                        frame = selection.query(f"  budget == {budget} & type == '{thing}' &  Qubits == {q} & g == {g}")
                        try:
                            desc = str(frame['description'].iloc[0])
                            sho = str(frame['shots'].iloc[0])
                            cop = str(frame['copies'].iloc[0])
                        except:
                            desc = 'no'
                            sho =  'no'
                            cop =  'no'
                        csv_writer.writerow([thing, budget, q, g,x,sho,cop,desc,result_type,f'Q:{q} g:{g}',q*q*g])
# %%
# tags = ['_LinGrow']
# title_dict = {'_LinGrow':'Linear Growth','_RNM_e1':'Scaled Noise (e_r=1)','_RNM_e0p1':'Scaled Noise (e_r=0.1)'}
# for tag in tags:
#     for QUBIT,g in zip([4,4,6,6,8,8,10],[1,16,1,16,1,16,1]):
#         fig, axes = plt.subplots(2,2,figsize=(15,10),sharey=False,sharex=True)
#         axes=axes.flatten()

#         sns.lineplot(data=t.query(
#             f'budget==0    & nlsp == 1 & Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0      &  type.str.endswith("{tag}")'),
#             ax=axes[0],x='copies',y='abs_error', hue='type',style='type' ,estimator='mean',err_style=None)
#         sns.lineplot(data=t.query(
#             f'budget==100000   & nlsp == 1 & Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0 &  type.str.endswith("{tag}")'),
#             ax=axes[1],x='copies',y='abs_error', hue='type',style='type' ,estimator='mean',err_style=None)
#         sns.lineplot(data=t.query(
#             f'budget==1000000  & nlsp == 1 & Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0 &  type.str.endswith("{tag}")'),
#             ax=axes[2],x='copies',y='abs_error', hue='type',style='type' ,estimator='mean',err_style=None)
#         sns.lineplot(data=t.query(
#             f'budget==10000000 & nlsp == 1 & Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0 &  type.str.endswith("{tag}")'),
#             ax=axes[3],x='copies',y='abs_error', hue='type',style='type' ,estimator='mean',err_style=None)

#         for i in axes:
#             i.set_yscale('log')
#             i.legend([],[] ,frameon=False)
#         axes[0].set_title('Infinite Budget')
#         axes[1].set_title('1M   Budget')
#         axes[2].set_title('10M  Budget')
#         axes[3].set_title('100M Budget')
#         handles, labels = axes[0].get_legend_handles_labels()
#         fig.legend(handles, labels, title='Mitigation Strategies:',loc='lower center', bbox_to_anchor=(0.5, -0.025),
#                 ncol=6, fancybox=True, shadow=True)
#         fig.suptitle(f"More shots for vn/CG at 10-9 {title_dict[tag]} Q={QUBIT} g={g}")
#         plt.tight_layout()
#         plt.savefig(f'./PLOTS/NO_CONST_{tag}_budgetsQ{QUBIT}g{g}.png')
#         plt.close(fig)
# %%
description = '3nlsp_half'
title_dict = {'_LinGrow':'Linear Growth'}#,'_RNM_e1':'Scaled Noise (e_r=1)','_RNM_e0p1':'Scaled Noise (e_r=0.1)'}
for  budget in [10**7,10**8,10**9,10**10]:
    fig, axes = plt.subplots(2,4,figsize=(15,10),sharey=False,sharex='col')
#    axes=axes.flatten():q
    plot_data=[]
    plot_data.append(scaled_panda.query(f'(description == "{description}" & budget=={budget}|budget==0)  & copies==1  & nlsp == 1  & abs_error != 0  &  (type.str.startswith("ZNE")|type.str.startswith("CDR")|type.str.startswith("vnCD"))'))
    plot_data.append(scaled_panda.query(f'(description == "{description}" & budget=={budget}|budget==0)  & copies==3  & nlsp == 1  & abs_error != 0  &  (type.str.startswith("VD")|type.str.startswith("vnCG")|type.str.startswith("CG"))'))

    plot_data=pd.concat(plot_data)
    for row, tag in enumerate(['_LinGrow']):
        for column,g in enumerate([1,16]):
                sns.lineplot(data=plot_data.query(
                f'    budget=={budget} & nlsp == 1 &  abs_error != 0  &  type.str.endswith("{tag}") & g == {g}'),
                ax = axes[column, row],x='Qubits',y='scaled_error', hue='type',style='type',estimator='mean',markers=marker_map,err_style=None,palette=cmap,
                     color={"VD_LinGrow": "pink", "vnCDR_LinGrow": "#742802"})
                axes[column, row].set_yscale('log')
                axes[column, row].legend([],[] ,frameon=False)
                axes[column, row].set_title(f'Model: {title_dict[tag]}, g = {g}')

    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, title='Mitigation Strategies:',loc='lower center', bbox_to_anchor=(0.5, -0.025),
            ncol=7, fancybox=True, shadow=True)
    if budget != 0:
        fig.suptitle(f"More shots CG/vnPerformance over qubits and depth at budget 10^{int(np.log10(budget))}  ")
        plt.savefig(f'./PLOTS/nlsp3_full-10^{int(np.log10(budget))}screens.png')

    else:
        fig.suptitle("Performance over qubits and depth at infinite shot limit.")
        plt.savefig('./PLOTS/nlsp3_full-budget_inf_shots_screens.png')
    plt.close(fig)


 # %%
base_folder = 'D:/Databases/VD_CDR/'
fig, axes = plt.subplots(6,2,figsize=(15,15),sharey=False,sharex='col')
#    axes=axes.flatten():q
description = 'hlaf_v_normal'
t['type'] = t['type'].str.replace(r'_LinGrow', '')

# plot_data=[]
# plot_data.append(scaled_panda.query(f'description == "{description}" & copies==1  & nlsp == 1  & abs_error != 0  &  type.str.startswith("ZNE")'))
# plot_data.append(scaled_panda.query(f'description == "{description}" &  copies==3  & nlsp == 1  & abs_error != 0  &  (type.str.startswith("VD")|type.str.startswith("vn")|type.str.startswith("CG"))'))
# plot_data=pd.concat(plot_data)

for row, budg in enumerate([10**5,10**6,10**7,10**8,10**9,10**10]):
    for column,g in enumerate([1,16]):
        plot_data=[]
        plot_data.append(t.query(f'budget == {budg} & description == "3nlsp_full" & copies==1  & nlsp == 1  & abs_error != 0  &  (type.str.startswith("ZNE"))'))
        plot_data.append(t.query(f'budget == {budg} & description == "3nlsp_full" & copies==3  & nlsp == 1  & abs_error != 0  &  (type.str.startswith("VD"))'))
        plot_data.append(t.query(f'budget == {budg*10} & description == "3nlsp_half" &  copies==3  & nlsp == 1  & abs_error != 0  &  type.str.startswith("vnCG")'))
        plot_data.append(t.query(f'budget == {budg*10} & description == "3nlsp_half" &  copies==1  & nlsp == 1  & abs_error != 0  &  type.str.startswith("vnCD")'))
        o = t.query(
    'abs_error != 0  & budget != 0 & description == "3nlsp_full" & copies == 1 & type.str.startswith("VD")').copy()
        o['type'] = 'Noisy'
        plot_data.append(o)
        plot_data=pd.concat(plot_data)

        sns.lineplot(data=plot_data.query(f' g == {g}'),
        ax = axes[row,column],x='Qubits',y='abs_error', hue='type',style='type',estimator='mean',markers=marker_map,err_style=None,palette=cmap,dashes=line_map)
        axes[row,column].set_yscale('log')
        axes[row,column].legend([],[] ,frameon=False)
        # axes[row,column].set_ylim([10**-2,10**1])
        # axes[row,1].set_yticklabels([])
        axes[row,1].set_ylabel('')
        axes[row,0].set_ylabel(r"$\overline{\langle \sigma_z \rangle}$ abs error")
        
        try:
            axes[row,column].set_title(f'Budget: 10^{int(np.log10(budg))}, g = {g}')
        except:
            axes[row,column].set_title(f'Budget: Infinite, g = {g}')



handles, labels = axes[-2][-2].get_legend_handles_labels()
fig.legend(handles, labels, title='Mitigation Strategies:',loc='lower center', bbox_to_anchor=(0.5, -0.025),
        ncol=7, fancybox=True, shadow=True)
if budg != 0:
    fig.suptitle(f" Comparison between budgets at different techniques")
    plt.savefig(base_folder+f'/PLOTS/A{description}_lineplot_budget_screens.png')

else:
    fig.suptitle("Comparison between budgets")
    plt.savefig(base_folder+f'/PLOTS/A{description}_lineplot_budget_screens.png')
plt.close(fig)
############ UNCOMMENT UP TO HERE

# %%
# %%
# tag = "_LinGrow"

# for QUBIT,g in zip([4,4,6,6,8,8,10],[1,16,1,16,1,16,10]):
#     fig, axes = plt.subplots(2,2,figsize=(15,10),sharey=False,sharex=True)
#     axes=axes.flatten()

#     sns.lineplot(data=t.query(
#         f'budget==0    & nlsp == 1 & Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0      &  type.str.endswith("{tag}")'),
#         ax=axes[0],x='copies',y='abs_error', hue='type' ,estimator='mean',err_style=None)
#     sns.lineplot(data=t.query(
#         f'budget==1000000   & nlsp == 1 & Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0 &  type.str.endswith("{tag}")'),
#         ax=axes[1],x='copies',y='abs_error', hue='type' ,estimator='mean',err_style=None)
#     sns.lineplot(data=t.query(
#         f'budget==10000000  & nlsp == 1 & Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0 &  type.str.endswith("{tag}")'),
#         ax=axes[2],x='copies',y='abs_error', hue='type' ,estimator='mean',err_style=None)
#     sns.lineplot(data=t.query(
#         f'budget==100000000 & nlsp == 1 & Qubits == {QUBIT} & depth == {QUBIT}*{g} & abs_error != 0 &  type.str.endswith("{tag}")'),
#         ax=axes[3],x='copies',y='abs_error', hue='type' ,estimator='mean',err_style=None)

#     for i in axes:
#         i.set_yscale('log')
#         i.legend([],[] ,frameon=False)
#     axes[0].set_title('Infinite Budget')
#     axes[1].set_title('1M   Budget')
#     axes[2].set_title('10M  Budget')
#     axes[3].set_title('100M Budget')
#     handles, labels = axes[0].get_legend_handles_labels()
#     fig.legend(handles, labels, title='Mitigation Strategies:',loc='lower center', bbox_to_anchor=(0.5, -0.025),
#             ncol=6, fancybox=True, shadow=True)
#     fig.suptitle(f"Scaled realistic noise with 0-filtered training data at {tag} Q={QUBIT} g={g}")
#     plt.tight_layout()
#     plt.savefig(f'./PLOTS/LINplot_0filtered_all_{tag}_budgetsQ{QUBIT}g{g}.png')
#     plt.savefig(f'./PLOTS/LINplot_0filtered_all_{tag}_budgetsQ{QUBIT}g{g}')
# # %% 
# ################# HIGH NOISE CIRCUITS CONSISTENT ##################################
# high_noise_abs, high_dnoise_rescaling, high_noise_values = load_multiple_files([4,4,6,6,8,8],[4,4*16,6,6*16,8,8*16],30,tags=['_wns'],shots=[None,100000],density_matrices=False,folder='./new_data/linearly_growing_noise/',nlsp_list = [1,2,3],train=True)
# # %%

# fig, axes = plt.subplots(3,2,figsize=(15,10),sharey=False,sharex=True)
# axes=axes.flatten()
# sns.lineplot(data=high_noise_abs.query('shots==0 & nlsp == 1 & Qubits == 4 & depth == 4 & ~type.str.startswith("CDR") & abs_error != 0'),x='copies',y='abs_error', hue='type' ,ax=axes[0 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_noise_abs.query('shots==0 & nlsp == 1 & Qubits == 6 & depth == 6 & ~type.str.startswith("CDR")& abs_error != 0 '),x='copies',y='abs_error', hue='type' ,ax=axes[2 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_noise_abs.query('shots==0 & nlsp == 1 & Qubits == 8 &  depth == 8 &~type.str.startswith("CDR") & abs_error != 0'),x='copies',y='abs_error', hue='type' ,ax=axes[4 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_noise_abs.query('shots==1 & nlsp == 1 & Qubits == 4 & depth == 4 & ~type.str.startswith("CDR") & abs_error != 0'),x='copies',y='abs_error', hue='type' ,ax=axes[1 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_noise_abs.query('shots==1 & nlsp == 1 & Qubits == 6 & depth == 6 & ~type.str.startswith("CDR")& abs_error != 0 '),x='copies',y='abs_error', hue='type' ,ax=axes[3 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_noise_abs.query('shots==1 & nlsp == 1 & Qubits == 8 &  depth == 8 &~type.str.startswith("CDR") & abs_error != 0'),x='copies',y='abs_error', hue='type' ,ax=axes[5 ],estimator='mean',err_style=None)

# axes[0].set_title('4 Qubits  inf shots')
# axes[2].set_title('6 Qubits  inf shots')
# axes[4].set_title('8 Qubits  inf shots')
# axes[1].set_title('4 Qubits  100K shots ')
# axes[3].set_title('6 Qubits  100K shots')
# axes[5].set_title('8 Qubits  100K shots')
# for i in axes:
#     i.set_yscale('log')
#     i.legend([],[] ,frameon=False)
# handles, labels = axes[0].get_legend_handles_labels()
# fig.legend(handles, labels, title='Mitigation Strategies:',loc='lower center', bbox_to_anchor=(0.5, -0.025),
#           ncol=4, fancybox=True, shadow=True)
# fig.suptitle("Low Noise (e_r=0.1) at g=1")
# plt.tight_layout()

# ################# PLOT  g=16
# fig, axes = plt.subplots(3,2,figsize=(15,10),sharey=False,sharex=True)
# axes=axes.flatten()
# sns.lineplot(data=high_noise_abs.query('shots==0 & nlsp == 1 & Qubits == 4 & depth == 4 *16& ~type.str.startswith("CDR") & ~type.str.startswith("Z")&~type.str.startswith("VD") & abs_error != 0'),x='copies',y='abs_error', hue='type' ,ax=axes[0 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_noise_abs.query('shots==0 & nlsp == 1 & Qubits == 6 & depth == 6 *16& ~type.str.startswith("CDR") & ~type.str.startswith("Z")&~type.str.startswith("VD")& abs_error != 0 '),x='copies',y='abs_error', hue='type' ,ax=axes[2 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_noise_abs.query('shots==0 & nlsp == 1 & Qubits == 8 &  depth == 8*16 &~type.str.startswith("CDR") & ~type.str.startswith("Z")&~type.str.startswith("VD") & abs_error != 0'),x='copies',y='abs_error', hue='type' ,ax=axes[4 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_noise_abs.query('shots==1 & nlsp == 1 & Qubits == 4 & depth == 4 *16& ~type.str.startswith("CDR") & ~type.str.startswith("Z")&~type.str.startswith("VD") & abs_error != 0'),x='copies',y='abs_error', hue='type' ,ax=axes[1 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_noise_abs.query('shots==1 & nlsp == 1 & Qubits == 6 & depth == 6 *16& ~type.str.startswith("CDR") & ~type.str.startswith("Z")&~type.str.startswith("VD")& abs_error != 0 '),x='copies',y='abs_error', hue='type' ,ax=axes[3 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_noise_abs.query('shots==1 & nlsp == 1 & Qubits == 8 &  depth == 8*16 &~type.str.startswith("CDR") & ~type.str.startswith("Z")&~type.str.startswith("VD") & abs_error != 0'),x='copies',y='abs_error', hue='type' ,ax=axes[5 ],estimator='mean',err_style=None)

# axes[0].set_title('4 Qubits  inf shots')
# axes[2].set_title('6 Qubits  inf shots')
# axes[4].set_title('8 Qubits  inf shots')
# axes[1].set_title('4 Qubits  100K shots ')
# axes[3].set_title('6 Qubits  100K shots')
# axes[5].set_title('8 Qubits  100K shots')
# for i in axes:
#     i.set_yscale('log')
#     i.legend([],[] ,frameon=False)
# handles, labels = axes[0].get_legend_handles_labels()
# fig.legend(handles, labels, title='Mitigation Strategies:',loc='lower center', bbox_to_anchor=(0.5, -0.025),
#           ncol=5, fancybox=True, shadow=True)
# fig.suptitle("Low Noise (e_r=0.1) at g=16")
# plt.tight_layout()
# # %%
# ####################################################### BOXPLOT #########################################################

# fig, axes = plt.subplots(3,2,figsize=(15,10),sharey=False,sharex=True)
# axes=axes.flatten()
# sns.boxplot(data=high_noise_abs.query('shots==0 & nlsp == 1 & Qubits == 4 & depth == 4 & ~type.str.startswith("CDR") & abs_error != 0 & copies == 6'),x='type',y='abs_error', hue='type' ,ax=axes[0 ])
# sns.boxplot(data=high_noise_abs.query('shots==0 & nlsp == 1 & Qubits == 6 & depth == 6 & ~type.str.startswith("CDR")& abs_error != 0  & copies == 6'),x='type',y='abs_error', hue='type' ,ax=axes[2 ])
# sns.boxplot(data=high_noise_abs.query('shots==0 & nlsp == 1 & Qubits == 8 &  depth == 8 &~type.str.startswith("CDR") & abs_error != 0 & copies == 6'),x='type',y='abs_error', hue='type' ,ax=axes[4 ])
# sns.boxplot(data=high_noise_abs.query('shots==1 & nlsp == 1 & Qubits == 4 & depth == 4 & ~type.str.startswith("CDR") & abs_error != 0 & copies == 6'),x='type',y='abs_error', hue='type' ,ax=axes[1 ])
# sns.boxplot(data=high_noise_abs.query('shots==1 & nlsp == 1 & Qubits == 6 & depth == 6 & ~type.str.startswith("CDR")& abs_error != 0  & copies == 6'),x='type',y='abs_error', hue='type' ,ax=axes[3 ])
# sns.boxplot(data=high_noise_abs.query('shots==1 & nlsp == 1 & Qubits == 8 &  depth == 8 &~type.str.startswith("CDR") & abs_error != 0 & copies == 6'),x='type',y='abs_error', hue='type' ,ax=axes[5 ])

# axes[0].set_title('4 Qubits  inf shots')
# axes[2].set_title('6 Qubits  inf shots')
# axes[4].set_title('8 Qubits  inf shots')
# axes[1].set_title('4 Qubits  100K shots ')
# axes[3].set_title('6 Qubits  100K shots')
# axes[5].set_title('8 Qubits  100K shots')
# for i in axes:
#     i.set_yscale('log')
#     i.legend([],[] ,frameon=False)
# handles, labels = axes[0].get_legend_handles_labels()
# fig.legend(handles, labels, title='Mitigation Strategies:',loc='lower center', bbox_to_anchor=(0.5, -0.025),
#           ncol=5, fancybox=True, shadow=True)
# fig.suptitle("Low Noise (e_r=0.1) at g=1")
# plt.tight_layout()

# ################# PLOT  g=16
# fig, axes = plt.subplots(3,2,figsize=(15,10),sharey=False,sharex=True)
# axes=axes.flatten()
# sns.boxplot(data=high_noise_abs.query('shots==0 & nlsp == 1 & Qubits == 4 & depth == 4 *16 & ~type.str.startswith("CDR") & abs_error != 0 & copies == 6'),x='type',y='abs_error', hue='type' ,ax=axes[0 ])
# sns.boxplot(data=high_noise_abs.query('shots==0 & nlsp == 1 & Qubits == 6 & depth == 6 *16 & ~type.str.startswith("CDR")& abs_error != 0  & copies == 6'),x='type',y='abs_error', hue='type' ,ax=axes[2 ])
# sns.boxplot(data=high_noise_abs.query('shots==0 & nlsp == 1 & Qubits == 8 &  depth == 8*16 & ~type.str.startswith("CDR") & abs_error != 0 & copies == 6'),x='type',y='abs_error', hue='type' ,ax=axes[4 ])
# sns.boxplot(data=high_noise_abs.query('shots==1 & nlsp == 1 & Qubits == 4 & depth == 4 *16 & ~type.str.startswith("CDR") & abs_error != 0 & copies == 6'),x='type',y='abs_error', hue='type' ,ax=axes[1 ])
# sns.boxplot(data=high_noise_abs.query('shots==1 & nlsp == 1 & Qubits == 6 & depth == 6 *16 & ~type.str.startswith("CDR")& abs_error != 0  & copies == 6'),x='type',y='abs_error', hue='type' ,ax=axes[3 ])
# sns.boxplot(data=high_noise_abs.query('shots==1 & nlsp == 1 & Qubits == 8 &  depth == 8*16 & ~type.str.startswith("CDR") & abs_error != 0 & copies == 6'),x='type',y='abs_error', hue='type' ,ax=axes[5 ])

# axes[0].set_title('4 Qubits  inf shots')
# axes[2].set_title('6 Qubits  inf shots')
# axes[4].set_title('8 Qubits  inf shots')
# axes[1].set_title('4 Qubits  100K shots ')
# axes[3].set_title('6 Qubits  100K shots')
# axes[5].set_title('8 Qubits  100K shots')
# for i in axes:
#     i.set_yscale('log')
#     i.legend([],[] ,frameon=False)
# handles, labels = axes[0].get_legend_handles_labels()
# fig.legend(handles, labels, title='Mitigation Strategies:',loc='lower center', bbox_to_anchor=(0.5, -0.025),
#           ncol=4, fancybox=True, shadow=True)
# fig.suptitle("Low Noise (e_r=0.1) at g=16")
# plt.tight_layout()
# # %%
# ################# HIGH DEPTH CIRCUITS CONSISTENT ##################################
# high_depth_abs, high_depth_rescaling, high_depth_values = load_multiple_files([4,6,8],[64,6*16,8*16],30,tags=['_RNM_e0-1','_SNM_e0-1'],shots=[None,100000],density_matrices=False,folder='./new_data/noise_model_comparison/',nlsp_list = [1,2,3],train=True)
# # %%
# fig, axes = plt.subplots(3,4,figsize=(15,10),sharey=False,sharex=True)
# axes=axes.flatten()
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & Qubits == 4 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[0 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & Qubits == 4 & type.str.endswith("SNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[1 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & Qubits == 4 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[2 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & Qubits == 4 & type.str.endswith("SNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[3 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & Qubits == 6 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[4 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & Qubits == 6 & type.str.endswith("SNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[5 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & Qubits == 6 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[6 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & Qubits == 6 & type.str.endswith("SNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[7 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & Qubits == 8 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[8 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & Qubits == 8 & type.str.endswith("SNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[9 ],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & Qubits == 8 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[10],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & Qubits == 8 & type.str.endswith("SNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[11],estimator='mean',err_style=None)


# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & Qubits == 4 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[0 ],estimator='mean',err_style=None,palette='hot')
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & Qubits == 4 & type.str.endswith("SNM_e0-1") & ~type.str.startswith("CDR") & type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[1 ],estimator='mean',err_style=None,palette='hot')
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & Qubits == 4 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[2 ],estimator='mean',err_style=None,palette='hot')
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & Qubits == 4 & type.str.endswith("SNM_e0-1") & ~type.str.startswith("CDR") & type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[3 ],estimator='mean',err_style=None,palette='hot')
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & Qubits == 6 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[4 ],estimator='mean',err_style=None,palette='hot')
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & Qubits == 6 & type.str.endswith("SNM_e0-1") & ~type.str.startswith("CDR") & type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[5 ],estimator='mean',err_style=None,palette='hot')
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & Qubits == 6 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[6 ],estimator='mean',err_style=None,palette='hot')
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & Qubits == 6 & type.str.endswith("SNM_e0-1") & ~type.str.startswith("CDR") & type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[7 ],estimator='mean',err_style=None,palette='hot')


# for i in axes:
#     i.set_yscale('log')
#     i.legend([],[] ,frameon=False)

# axes[0].set_title('4 Qubits RNM inf. ')
# axes[1].set_title('4 Qubits SNM inf. ')
# axes[2].set_title('4 Qubits RNM 100K ')
# axes[3].set_title('4 Qubits SNM 100K ')
# axes[4].set_title('6 Qubits RNM inf  ')
# axes[5].set_title('6 Qubits SNM inf  ')
# axes[6].set_title('6 Qubits RNM 100k  ')
# axes[7].set_title('6 Qubits SNM 100k  ')
# axes[8 ].set_title('8 Qubits RNM inf  ')
# axes[9 ].set_title('8 Qubits SNM inf  ')
# axes[10].set_title('8 Qubits RNM 100k  ')
# axes[11].set_title('8 Qubits SNM 100k  ')
# handles, labels = axes[0].get_legend_handles_labels()
# fig.legend(handles, ['CGVD', 'VD', 'vnCDR', 'vnCGVD'], title='Mitigation Strategies:',loc='lower center', bbox_to_anchor=(0.5, -0.025),
#           ncol=4, fancybox=True, shadow=True)
# fig.suptitle("Noise model comparison at g=16")
# plt.tight_layout()
# # %%
# fig, axes = plt.subplots(3,4,figsize=(15,10),sharey=False,sharex="col")
# axes=axes.flatten()

# VD_df = high_depth_abs.query(' nlsp == 1 & copies == 3 & type.str.startswith("VD")')
# nVD_df = high_depth_abs.query(' nlsp == 1 & copies == 6 & ~type.str.startswith("VD")')
# plot_stuff = pd.concat([VD_df,nVD_df])
# sns.boxplot(data=plot_stuff.query('shots==0 & nlsp == 1 & Qubits == 4 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[0 ],showfliers = False)
# sns.boxplot(data=plot_stuff.query('shots==0 & nlsp == 1 & Qubits == 4 & type.str.endswith("SNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[1 ],showfliers = False)
# sns.boxplot(data=plot_stuff.query('shots==1 & nlsp == 1 & Qubits == 4 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[2 ],showfliers = False)
# sns.boxplot(data=plot_stuff.query('shots==1 & nlsp == 1 & Qubits == 4 & type.str.endswith("SNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[3 ],showfliers = False)
# sns.boxplot(data=plot_stuff.query('shots==0 & nlsp == 1 & Qubits == 6 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[4 ],showfliers = False)
# sns.boxplot(data=plot_stuff.query('shots==0 & nlsp == 1 & Qubits == 6 & type.str.endswith("SNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[5 ],showfliers = False)
# sns.boxplot(data=plot_stuff.query('shots==1 & nlsp == 1 & Qubits == 6 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[6 ],showfliers = False)
# sns.boxplot(data=plot_stuff.query('shots==1 & nlsp == 1 & Qubits == 6 & type.str.endswith("SNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[7 ],showfliers = False)
# sns.boxplot(data=plot_stuff.query('shots==0 & nlsp == 1 & Qubits == 8 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[8 ],showfliers = False)
# sns.boxplot(data=plot_stuff.query('shots==0 & nlsp == 1 & Qubits == 8 & type.str.endswith("SNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[9 ],showfliers = False)
# sns.boxplot(data=plot_stuff.query('shots==1 & nlsp == 1 & Qubits == 8 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[10],showfliers = False)
# sns.boxplot(data=plot_stuff.query('shots==1 & nlsp == 1 & Qubits == 8 & type.str.endswith("SNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[11],showfliers = False)


# for i in axes:
#     i.set_yscale('log')
#     i.legend([],[] ,frameon=False)
#     i.set_yticks(range(4)) # <--- set the ticks first

#     i.set_yticklabels(['','','',''])

# axes[0].set_title('4 Qubits RNM inf. ')
# axes[1].set_title('4 Qubits SNM inf. ')
# axes[2].set_title('4 Qubits RNM 100K ')
# axes[3].set_title('4 Qubits SNM 100K ')
# axes[4].set_title('6 Qubits RNM inf  ')
# axes[5].set_title('6 Qubits SNM inf  ')
# axes[6].set_title('6 Qubits RNM 100k  ')
# axes[7].set_title('6 Qubits SNM 100k  ')
# axes[8 ].set_title('8 Qubits RNM inf  ')
# axes[9 ].set_title('8 Qubits SNM inf  ')
# axes[10].set_title('8 Qubits RNM 100k  ')
# axes[11].set_title('8 Qubits SNM 100k  ')
# handles, labels = axes[0].get_legend_handles_labels()
# fig.legend(handles, ['VD','CGVD','vnCDR','vnCGVD'], title='Mitigation Strategies:',loc='lower center', bbox_to_anchor=(0.5, -0.025),
#           ncol=4, fancybox=True, shadow=True)
# fig.suptitle("Noise model comparison at g=16")
# plt.tight_layout()
# # %%
# ######### PLOT WITH ONLY RNM
# fig, axes = plt.subplots(3,2,figsize=(10,10),sharey="col",sharex="col")
# axes=axes.flatten()

# VD_df = high_depth_abs.query(' nlsp == 1 & copies == 3 & type.str.startswith("VD")')
# nVD_df = high_depth_abs.query(' nlsp == 1 & copies == 6 & ~type.str.startswith("VD")')
# plot_stuff = pd.concat([VD_df,nVD_df])
# sns.boxplot(data=plot_stuff.query('shots==0 & nlsp == 1 & Qubits == 4 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[0 //2],showfliers = False)
# sns.boxplot(data=plot_stuff.query('shots==1 & nlsp == 1 & Qubits == 4 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[2 //2],showfliers = False)
# sns.boxplot(data=plot_stuff.query('shots==0 & nlsp == 1 & Qubits == 6 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[4 //2],showfliers = False)
# sns.boxplot(data=plot_stuff.query('shots==1 & nlsp == 1 & Qubits == 6 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[6 //2],showfliers = False)
# sns.boxplot(data=plot_stuff.query('shots==0 & nlsp == 1 & Qubits == 8 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[8 //2],showfliers = False)
# sns.boxplot(data=plot_stuff.query('shots==1 & nlsp == 1 & Qubits == 8 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR")'),x='type',y='abs_error',hue='type' ,ax=axes[10//2],showfliers = False)


# for i in axes:
#     i.set_yscale('log')
#     i.set_xlabel('')
#     i.legend([],[] ,frameon=False)
#     i.set_xticks(range(4)) # <--- set the ticks first
#     i.set_xticklabels(['','','',''])
#     #i.set_xticklabels(['VD','CGVD','vnCDR','vnCGVD'])

# axes[0 //2].set_title('4 Qubits RNM inf. ')
# axes[2 //2].set_title('4 Qubits RNM 100K ')
# axes[4 //2].set_title('6 Qubits RNM inf  ')
# axes[6 //2].set_title('6 Qubits RNM 100k  ')
# axes[8 //2].set_title('8 Qubits RNM inf  ')
# axes[10//2].set_title('8 Qubits RNM 100k  ')
# handles, labels = axes[0].get_legend_handles_labels()
# fig.legend(handles, ['VD','CGVD','vnCDR','vnCGVD'], title='Mitigation Strategies:',loc='lower center', bbox_to_anchor=(0.5, -0.04),
#           ncol=4, fancybox=True, shadow=True)
# fig.suptitle("Realistic noise at e_r=0.1, g=16, optimal copies")
# plt.tight_layout()
#  # %%
#  # %%
# fig, axes = plt.subplots(3,2,figsize=(10,10),sharey='col',sharex=True)
# axes=axes.flatten()
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & Qubits == 4 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[0 //2],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & Qubits == 4 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[2 //2],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & Qubits == 6 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[4 //2],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & Qubits == 6 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[6 //2],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & Qubits == 8 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[8 //2],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & Qubits == 8 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[10//2],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & Qubits == 4 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & type.str.startswith("vn")' ),x='copies',y='abs_error', hue='type' , ax=axes[0 //2],estimator='mean',err_style=None,palette='hot')
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & Qubits == 4 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & type.str.startswith("vn")' ),x='copies',y='abs_error', hue='type' , ax=axes[2 //2],estimator='mean',err_style=None,palette='hot')
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & Qubits == 6 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & type.str.startswith("vn")' ),x='copies',y='abs_error', hue='type' , ax=axes[4 //2],estimator='mean',err_style=None,palette='hot')
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & Qubits == 6 & type.str.endswith("RNM_e0-1") & ~type.str.startswith("CDR") & type.str.startswith("vn")' ),x='copies',y='abs_error', hue='type' , ax=axes[6 //2],estimator='mean',err_style=None,palette='hot')


# for i in axes:
#     i.set_yscale('log')
#     i.legend([],[], frameon=False)

# axes[0 //2].set_title('4 Qubits RNM inf. ')
# axes[2 //2].set_title('4 Qubits RNM 100K ')
# axes[4 //2].set_title('6 Qubits RNM inf  ')
# axes[6 //2].set_title('6 Qubits RNM 100k  ')
# axes[8 //2].set_title('8 Qubits RNM inf  ')
# axes[10//2].set_title('8 Qubits RNM 100k  ')
# handles, labels = axes[0].get_legend_handles_labels()
# fig.legend(handles, ['CGVD', 'VD', 'vnCDR', 'vnCGVD'], title='Mitigation Strategies:',loc='lower center', bbox_to_anchor=(0.5, -0.025),
#           ncol=4, fancybox=True, shadow=True)
# fig.suptitle("g=16 Performance")
# plt.tight_layout()

#  # %%
 
# ######################## HIGH DEPTH CIRCUITS ####################################
# high_depth_abs, high_depth_rescaling, high_depth_values = load_multiple_files([4,4],[4,64],30,tags=['_RNM_e0-128'],shots=[None,100000],density_matrices=False,folder='./new_data/high_depth_runs/',nlsp_list = [1,2,3],train=True)
# ####################### SImple Noise model ###################################
# # %%


# fig, axes = plt.subplots(2,2,figsize=(10,10),sharey="row",sharex=True)
# axes=axes.flatten()
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & depth == 4   & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[0//2],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & depth == 4   & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[4//2],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & depth == 64  & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[2//2],estimator='mean',err_style=None)
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & depth == 64  & ~type.str.startswith("CDR") & ~type.str.startswith("vn")'),x='copies',y='abs_error', hue='type' ,ax=axes[6//2],estimator='mean',err_style=None)

# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & depth == 4   & ~type.str.startswith("CDR") & type.str.startswith("vn")' ),x='copies',y='abs_error', hue='type' ,ax=axes[0//2],estimator='mean',err_style=None,palette='hot')
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & depth == 4   & ~type.str.startswith("CDR") & type.str.startswith("vn")' ),x='copies',y='abs_error', hue='type' ,ax=axes[4//2],estimator='mean',err_style=None,palette='hot')
# sns.lineplot(data=high_depth_abs.query('shots==0 & nlsp == 1 & depth == 64  & ~type.str.startswith("CDR") & type.str.startswith("vn")' ),x='copies',y='abs_error', hue='type' ,ax=axes[2//2],estimator='mean',err_style=None,palette='hot')
# sns.lineplot(data=high_depth_abs.query('shots==1 & nlsp == 1 & depth == 64  & ~type.str.startswith("CDR") & type.str.startswith("vn")' ),x='copies',y='abs_error', hue='type' ,ax=axes[6//2],estimator='mean',err_style=None,palette='hot')


# for i in axes:
#     i.set_yscale('log')
#     i.legend([],[], frameon=False)

# axes[0//2].set_title('4 Qubits RNM inf. shots - depth 4')
# axes[4//2].set_title('4 Qubits RNM 100K shots - depth 4')
# axes[2//2].set_title('4 Qubits RNM at inf shots - depth 64')
# axes[6//2].set_title('4 Qubits RNM at 100k shots - depth 64')
# handles, labels = axes[0].get_legend_handles_labels()
# fig.legend(handles, ['CGVD', 'VD', 'vnCDR', 'vnCGVD'], title='Mitigation Strategies:',loc='lower center', bbox_to_anchor=(0.5, -0.025),
#           ncol=4, fancybox=True, shadow=True)
# fig.suptitle("Performance at varying depth")
# plt.tight_layout()######################## VARIABLE DEPTH BIG CIRCUITS ####################################
# # %%
# qubits = [4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10]
# depths = [q*[1,4,16][i%3] for i,q in enumerate(qubits) ]
# test_abs, test_rescaling, values = load_multiple_files(qubits[0:-3],depths,30,tags=['_RNM'],shots=[None],density_matrices=True,folder='./new_data/VD_4-10_scaling/',nlsp_list = [1],train=False)
# # %%
# simple_abs, simple_rescaling, simple_values = load_multiple_files(qubits[0:-3],depths,30,tags=['_SNM'],shots=[None],density_matrices=False,folder='./new_data/VD_4-10_scaling/',nlsp_list = [1],train=False)
# # %%
# test_abs.Qubits = test_abs.Qubits.astype(str)
# values['g']=values['depth']//values['Qubits'].astype(int)

# test_rescaling['g']=test_rescaling['depth']//test_rescaling['Qubits'].astype(int)
# fig, axes = plt.subplots(2,2,figsize=(10,10),sharey=False,sharex=True)
# axes = axes.flatten()
# #sns.boxplot(data=test_abs.query('type.str.startswith("coheren")  & copies == 1'),hue='g',x='Qubits',y='abs_error',ax=axes,palette='Set2')
# sns.lineplot(data=test_abs.query('type.str.startswith("coheren")  & copies == 1'),hue='g',x='Qubits',y='abs_error',ax=axes[0],estimator='mean',err_style=None,marker="o",palette='Set2')
# sns.scatterplot(data=test_abs.query('type.str.startswith("coheren")  & copies == 1'),hue='g',x='Qubits',y='abs_error',ax=axes[0],alpha=0.2,palette='Set2')
# sns.lineplot(data=test_abs.query('type.str.startswith(   "noise")  & copies == 2'),hue='g',x='Qubits',y='abs_error',ax=axes[1],estimator='mean',err_style=None,marker="o",palette='Set2')
# # sns.scatterplot(data=test_abs.query('type.str.startswith("noise")  & copies == 2'),hue='g',x='Qubits',y='abs_error',ax=axes[1],alpha=0.2,palette='Set2')
# sns.lineplot(data=test_abs.query('type.str.startswith(   "eigen")  & copies == 2'),hue='g',x='Qubits',y='abs_error',ax=axes[2],estimator='mean',err_style=None,marker="o",palette='Set2')
# sns.scatterplot(data=test_abs.query('type.str.startswith("eigen")  & copies == 3'),hue='g',x='Qubits',y='abs_error',ax=axes[2],alpha=0.2,palette='Set2')
# sns.lineplot(data=test_abs.query('type.str.startswith(   "VD")  & copies ==    3'),hue='g',x='Qubits',y='abs_error',ax=axes[3],estimator='mean',err_style=None,marker="o",palette='Set2')
# # sns.scatterplot(data=test_abs.query('type.str.startswith("VD")  & copies == 1'),hue='g',x='Qubits',y='abs_error',ax=axes[3],alpha=0.2,palette='Set2')

# for i in axes:
#     i.set_yscale('log')
#     i.legend([],[], frameon=False)

# axes[0].set_title('Coherence Mismatch')
# axes[1].set_title('Noise Floor')
# axes[2].set_title('Eigenvalue Ratio')
# axes[3].set_title('Virtual Distillation (optimal result)')
# axes[2].set_ylabel('Ratio of subdominant/dominant noisy eigenvectors')

# handles, labels = axes[0].get_legend_handles_labels()
# fig.legend(handles, labels, title='depth scaling (g):',loc='lower center', bbox_to_anchor=(0.5, -0.035),
#           ncol=6, fancybox=True, shadow=True)
# fig.suptitle("Analytical quantities")
# plt.tight_layout()

# # %%
# # test_rescaling_2c= test_rescaling.query('copies == 2')
# # simple_rescaling_2c= simple_rescaling.query('copies == 2')
# test_rescaling_3c= simple_rescaling.query('copies == 2')
# # simple_rescaling_3c= simple_rescaling.query('copies == 3')

# # test_rescaling_2c['g']=test_rescaling_2c['depth']//test_rescaling_2c['Qubits']
# # simple_rescaling_2c['g']=simple_rescaling_2c['depth']//simple_rescaling_2c['Qubits']
# test_rescaling_3c['g']=test_rescaling_3c['depth']//test_rescaling_3c['Qubits']
# # simple_rescaling_3c['g']=simple_rescaling_3c['depth']//simple_rescaling_3c['Qubits']
# # %%

# fig, axes = plt.subplots(1,1,figsize=(10,10),sharey=True)
# sns.lineplot(data=test_rescaling_3c,x='Qubits',y='error rescaling factor',hue='g',ax = axes,marker='o')
# axes.set_title('Simple noise model; error rate = 0.5')
# axes.set_yscale('log')



# # %%
# #%%
# fig, axes = plt.subplots(2,2,figsize=(15,10),sharey=True)
# axes=axes.flatten()


# # sns.lineplot(data=test_rescaling_q,x='Qubits',y='error rescaling factor',hue='depth',ax=axes[1-1],estimator='mean',err_style=None)
# # sns.lineplot(data=simple_rescaling_q,x='Qubits',y='error rescaling factor',hue='depth',ax=axes[2-1],estimator='mean',err_style=None)

# # 3 copies
# simple_rescaling_g1 = [0.068839,0.035140,0.024319]
# simple_rescaling_g2 = [0.020221,0.014431,0.005401]
# simple_rescaling_g3 = [0.011885,0.004727,0.002239]
# # 2 copies
# normal_1 = [0.105377,0.186793,0.119646]
# normal_2 = [0.104954,0.097419,0.062211]
# normal_3  = [0.043630,0.011829,0.008805]

# normal3_1 = [0.122619,0.194748,	0.124864]
# normal3_2 = [0.110335,0.099266,0.063639]
# normal3_3  = [0.043794	,0.011491	,0.009144]

# simple_rescaling2_g1 = [0.070403,0.028995	,0.022614]
# simple_rescaling2_g2 = [0.063509,0.019331,0.005461]
# simple_rescaling2_g3 = [0.059970	,0.015798,0.004446]
# y = [4,6,8]

# axes[0].plot(y,simple_rescaling2_g1,label='G=1', linestyle='--', marker='o')
# axes[0].plot(y,simple_rescaling2_g2,label='G=4', linestyle='--', marker='o')
# axes[0].plot(y,simple_rescaling2_g3,label='G=6', linestyle='--', marker='o')
# axes[3].set_xticks([4,6,8])

# axes[2].plot(y,normal_1,label='G=1', linestyle='--', marker='o')
# axes[2].plot(y,normal_2,label='G=4', linestyle='--', marker='o')
# axes[2].plot(y,normal_3,label='G=6', linestyle='--', marker='o')
# axes[2].set_xticks([4,6,8])

# axes[0+1].plot(y,simple_rescaling_g1,label='G=1', linestyle='--', marker='o')
# axes[0+1].plot(y,simple_rescaling_g2,label='G=4', linestyle='--', marker='o')
# axes[0+1].plot(y,simple_rescaling_g3,label='G=6', linestyle='--', marker='o')
# axes[1].set_xticks([4,6,8])

# axes[2+1].plot(y,normal3_1,label='G=1', linestyle='--', marker='o')
# axes[2+1].plot(y,normal3_2,label='G=4', linestyle='--', marker='o')
# axes[2+1].plot(y,normal3_3,label='G=6', linestyle='--', marker='o')
# for i in axes:
#     i.set_yscale('log')
# axes[0].set_title('Simple - 2 copies')
# axes[1].set_title('Simple - 3 copies')
# axes[2].set_title('Normal - 2 copies')
# axes[3].set_title('normal - 3 copies')
# axes[0].set_ylabel('error rescaling factor')
# axes[2].set_ylabel('error rescaling factor')
# axes[2].set_xlabel('Qubits')
# axes[3].set_xlabel('Qubits')
# axes[0].set_xticks([4,6,8])

# axes[0].legend()

# axes[1].legend()


# # %%


# abs_q = test_abs.query('nlsp == 1 & shots == 0')
# abs_simp = simple_abs.query('nlsp == 1 & shots == 0')

# fig, axes = plt.subplots(3,2,figsize=(10,10),sharey=True)
# axes=axes.flatten()
# sns.boxplot(data=simple_abs.query('Qubits==4 & type.str.endswith("wns") & ~type.str.startswith("noise") & ~type.str.startswith("coherent") & ~type.str.startswith("eigen")'),x='depth',y='abs_error',hue='copies',ax=axes[0])
# sns.boxplot(data=simple_abs.query('Qubits==6 & type.str.endswith("wns") & ~type.str.startswith("noise") & ~type.str.startswith("coherent") & ~type.str.startswith("eigen")'),x='depth',y='abs_error',hue='copies',ax=axes[2])
# sns.boxplot(data=simple_abs.query('Qubits==8 & type.str.endswith("wns") & ~type.str.startswith("noise") & ~type.str.startswith("coherent") & ~type.str.startswith("eigen")'),x='depth',y='abs_error',hue='copies',ax=axes[4])

# sns.boxplot(data=test_abs.query('Qubits==4 & type.str.endswith("wns") & ~type.str.startswith("noise") & ~type.str.startswith("coherent") & ~type.str.startswith("eigen")'),x='depth',y='abs_error',hue='copies',ax=axes[1])
# sns.boxplot(data=test_abs.query('Qubits==6 & type.str.endswith("wns") & ~type.str.startswith("noise") & ~type.str.startswith("coherent") & ~type.str.startswith("eigen")'),x='depth',y='abs_error',hue='copies',ax=axes[3])
# sns.boxplot(data=test_abs.query('Qubits==8 & type.str.endswith("wns") & ~type.str.startswith("noise") & ~type.str.startswith("coherent") & ~type.str.startswith("eigen")'),x='depth',y='abs_error',hue='copies',ax=axes[5])


# axes[0].set_title('qubit 4,  normal noise model WNS inf shots')
# axes[2].set_title('qubit 6,  normal noise model WNS inf shots')
# axes[4].set_title('qubit 8,  normal noise model WNS inf shots')
# axes[1].set_title('qubit 4,  simple noise model WNS inf shots')
# axes[3].set_title('qubit 6,  simple noise model WNS inf shots')
# axes[5].set_title('qubit 8,  simple noise model WNS inf shots')

# for i in axes:
#     i.set_yscale('log')
#     i.legend([],[], frameon=False)
# axes[-1].legend(frameon=True,bbox_to_anchor=(1.02, 0.15), loc='upper left', borderaxespad=0)
# plt.tight_layout()
# # %% 
# fig, axes = plt.subplots(3,2,figsize=(10,10),sharey=True)

# axes=axes.flatten()

# for i,q in enumerate([4,5,6,7,8,9]):
#     sns.boxplot(data=simple_abs.query('Qubits=={} & (type.str.startswith("coherent")|type.str.startswith("noise"))'.format(q)),x='depth',y='abs_error',hue='type',ax=axes[0+2*i])
#     sns.boxplot(data=test_abs.query('Qubits=={} & (type.str.startswith("coherent")|type.str.startswith("noise"))'.format(q)),x='depth',y='abs_error',hue='type',ax=axes[1+2*i])
#     axes[2*i].set_title('{} qubits, NNM '.format(q))
#     axes[1+2*i].set_title('{} qubits, SNM '.format(q))
    
# for i in axes:
#     i.set_yscale('log')
#     i.legend([],[], frameon=False)
# axes[-1].legend(frameon=True,bbox_to_anchor=(1.02, 0.15), loc='upper left', borderaxespad=0)
# plt.tight_layout()

# # %%
# # %%


# # sns.boxplot(data=simple_abs.query('Qubits==4 & (type.str.startswith("coherent")|type.str.startswith("noise"))'),x='depth',y='abs_error',hue='type',ax=axes[0])
# # sns.boxplot(data=simple_abs.query('Qubits==6 & (type.str.startswith("coherent")|type.str.startswith("noise"))'),x='depth',y='abs_error',hue='type',ax=axes[2])
# # sns.boxplot(data=simple_abs.query('Qubits==8 & (type.str.startswith("coherent")|type.str.startswith("noise"))'),x='depth',y='abs_error',hue='type',ax=axes[4])
# # sns.boxplot(data=test_abs.query('Qubits==4 & (type.str.startswith("coherent")|type.str.startswith("noise"))'),x='depth',y='abs_error',hue='type',ax=axes[1])
# # sns.boxplot(data=test_abs.query('Qubits==6 & (type.str.startswith("coherent")|type.str.startswith("noise"))'),x='depth',y='abs_error',hue='type',ax=axes[3])
# # sns.boxplot(data=test_abs.query('Qubits==8 & (type.str.startswith("coherent")|type.str.startswith("noise"))'),x='depth',y='abs_error',hue='type',ax=axes[5])

# axes[0].set_title('4 qubits, NNM Coherence mismatch & noise floor')

# for i in axes:
#     i.set_yscale('log')
#     i.legend([],[], frameon=False)
# axes[-1].legend(frameon=True,bbox_to_anchor=(1.02, 0.15), loc='upper left', borderaxespad=0)
# plt.tight_layout()
# #%%

# #%%
# # VD_CDR(8,8,30,10,100,6,[1,2,3],shots=[None],
#                             #    folder='./new_data/100_factor_rescaled/',autoprocessing=True,plots=True,only_stat_plots=True,
#                             #    density_matrices = False,extra_tags=['_wns'])# %%
# ############# VARIABLE DEPTH CIRCUITS ########################
# test_abs, test_rescaling, values = load_multiple_files([4,4,4,4,6,6,6,6],[4,6,8,10,4,6,8,10],10,tags=['','_wns'],shots=[None,100000],density_matrices=False,folder='./new_data/variable_depth/')
# # %%


# abs_q = test_abs.query('nlsp == 1 & shots == 0 & copies == 6')
# fig, axes = plt.subplots(2,2,figsize=(10,10),sharey=True)
# axes=axes.flatten()
# sns.boxplot(data=abs_q.query('Qubits==4 & ~type.str.endswith("wns") & ~type.str.startswith("noise") & ~type.str.startswith("coherent") & ~type.str.startswith("eigen")'),x='depth',y='abs_error',hue='type',ax=axes[0])
# sns.boxplot(data=abs_q.query('Qubits==4 &  type.str.endswith("wns") & ~type.str.startswith("noise") & ~type.str.startswith("coherent") & ~type.str.startswith("eigen")'),x='depth',y='abs_error',hue='type',ax=axes[1])
# sns.boxplot(data=abs_q.query('Qubits==6 & ~type.str.endswith("wns") & ~type.str.startswith("noise") & ~type.str.startswith("coherent") & ~type.str.startswith("eigen")'),x='depth',y='abs_error',hue='type',ax=axes[2])
# sns.boxplot(data=abs_q.query('Qubits==6 &  type.str.endswith("wns") & ~type.str.startswith("noise") & ~type.str.startswith("coherent") & ~type.str.startswith("eigen")'),x='depth',y='abs_error',hue='type',ax=axes[3])

# axes[0].set_title('qubit 4, 6 copies , normal noise + 100K shots')
# axes[1].set_title('qubit 4, 6 copies , wns + 100K shots')
# axes[2].set_title('qubit 6, 6 copies , normal noise + 100K shots')
# axes[3].set_title('qubit 6, 6 copies , wns + 100K shots')

# for i in axes:
#     i.set_yscale('log')
#     i.legend([],[], frameon=False)
# axes[-1].legend(frameon=True,bbox_to_anchor=(1.02, 0.15), loc='upper left', borderaxespad=0)
# # %%



# # %%
# #################### COHERENT MISMATCH 
# abs_q = test_abs.query('nlsp == 1 & shots == 0 & copies == 2')
# fig, axes = plt.subplots(1,2,figsize=(10,5),sharey=True)
# axes=axes.flatten()
# sns.boxplot(data=abs_q.query('Qubits==4 & (type.str.startswith("coherent")|type.str.startswith("noise"))'),x='depth',y='abs_error',hue='type',ax=axes[0])
# sns.boxplot(data=abs_q.query('Qubits==6 & (type.str.startswith("coherent")|type.str.startswith("noise"))'),x='depth',y='abs_error',hue='type',ax=axes[1])

# axes[0].set_title('4 qubits, coherence mismatch vs depth')
# axes[1].set_title('6 qubits, coherence mismatch vs depth')

# for i in axes:
#     i.set_yscale('log')
#     i.legend([],[], frameon=False)
# axes[-1].legend(frameon=True,bbox_to_anchor=(1.02, 0.15), loc='upper left', borderaxespad=0)


# # %%
# ##################### Eigenenergy ratio


# abs_q = test_abs.query('nlsp == 1 & shots == 0 & copies == 3')
# fig, axes = plt.subplots(2,1,figsize=(5,7),sharey=False)
# axes=axes.flatten()
# sns.boxplot(data=abs_q.query(' ~type.str.endswith("wns") & type.str.startswith("eigen")'),x='depth',y='abs_error',hue='Qubits',ax=axes[0])
# sns.boxplot(data=abs_q.query('  type.str.endswith("wns") & type.str.startswith("eigen")'),x='depth',y='abs_error',hue='Qubits',ax=axes[1])

# axes[0].set_title('Subominant noisy eigenenergy over dominant - without noise scaling')
# axes[1].set_title('With noise scaling')

# for i in axes:
#     i.set_yscale('log')
#     i.legend([],[], frameon=False)
# axes[-1].legend(frameon=True,bbox_to_anchor=(1.02, 0.15), loc='upper left', borderaxespad=0)

# # density_matrix_test_abs, density_matrix_test_rescaling,dfs_standard_values = load_multiple_files([8],[8],30,tags=['_wns'],density_matrices=False,folder='./new_data/thirty_seeds/')
# # %%
# density_matrix_test_abs, density_matrix_test_rescaling,dfs_standard_values = load_multiple_files([4,6,8],[4,6,8],30,tags=['_wns'],density_matrices=False,folder='./new_data/100_factor_rescaled/')
# # %%
# dfs= dfs_standard_values
# data_floor=dfs.loc[(((dfs['type']== 'noise_floor_wns') & (dfs["nlsp"]==1) & (dfs["shots"]==0))& (dfs["Qubits"]==4))]
# data_VD=dfs.loc[(((dfs['type']== 'VD_wns') & (dfs["nlsp"]==1) & (dfs["shots"]==0))& (dfs["Qubits"]==4))]

# mean_VD=data_VD.groupby('copies').mean()
# mean_floor= data_floor.groupby('copies').mean()
# figure_5 = np.abs(mean_floor-mean_VD)

# fig, axes = plt.subplots(1,1,figsize=(10,10))
# sns.lineplot(data=figure_5,ax=axes)
# axes.set_yscale('log')
# # %%
# density_matrix_test_abs, density_matrix_test_rescaling,dfs_standard_values = load_multiple_files([4,6,8],[5],10,tags=['','_wns'],density_matrices=False,folder='./new_data/twenty_seeds/')
# # %%
# dfs= dfs_standard_values
# data_floor=dfs.loc[(((dfs['type']== 'noise_floor_wns') & (dfs["nlsp"]==1) & (dfs["shots"]==0))& (dfs["Qubits"]==4))]
# data_VD=dfs.loc[(((dfs['type']== 'VD_wns') & (dfs["nlsp"]==1) & (dfs["shots"]==0))& (dfs["Qubits"]==4))]

# mean_VD=data_VD.groupby('copies').mean()
# mean_floor= data_floor.groupby('copies').mean()
# figure_5 = np.abs(mean_floor-mean_VD)

# fig, axes = plt.subplots(1,1,figsize=(10,10))
# sns.lineplot(data=figure_5,ax=axes)
# axes.set_yscale('log')
# # %%
# density_matrix_test_abs_high, density_matrix_test_rescaling_high,dfs_standard_values_high = load_multiple_files([4,6,8,10],[4,6,8,10],30,tags=['_wns'],density_matrices=False,folder='./new_data/100_factor_rescaled/',shots=[None,10000,100000,1000000])
# density_matrix_test_abs_low, density_matrix_test_rescaling_low,dfs_standard_values_low = load_multiple_files([4,6,8,10],[4,6,8,10],30,tags=['_wns'],density_matrices=False,folder='./new_data/thirty_seeds/',shots=[None,100000])
# # %%

# df_plot_high = density_matrix_test_abs_high.query('nlsp == 1 & shots == 0 & copies == 6')
# df_plot_high3 = density_matrix_test_abs_high.query('nlsp == 1 & shots == 1 & copies == 6')
# df_plot_low  = density_matrix_test_abs_high.query('nlsp == 1 & shots == 2 & copies == 6')
# df_plot_low3 = density_matrix_test_abs_high.query('nlsp == 1 & shots == 3 & copies == 6')

# fig, axes = plt.subplots(2,2,figsize=(15,10),sharey=True)
# axes=axes.flatten()
# sns.boxplot(data=df_plot_high,x='type',y='abs_error',hue='Qubits',ax=axes[0])
# sns.boxplot(data=df_plot_low,x='type',y='abs_error',hue='Qubits',ax=axes[2])
# sns.boxplot(data=df_plot_high3,x='type',y='abs_error',hue='Qubits',ax=axes[1])
# sns.boxplot(data=df_plot_low3,x='type',y='abs_error',hue='Qubits',ax=axes[3])



# axes[0].set_title('HIGH noise infinite shots - 6 copies, nlsp 1')
# axes[0].set_yscale('log')
# axes[2].set_title('HIGH noise 10k shots - 6 copies, nlsp 1')
# axes[2].set_yscale('log')
# axes[1].set_title(' HIGH noise 100K shots - 6 copies, nlsp 1')
# axes[1].set_yscale('log')
# axes[3].set_title(' HIGH noise 1000k shots - 6 copies, nlsp 1')
# axes[3].set_yscale('log')

# # %%
# df_plot_0 = density_matrix_test_abs_high.query('nlsp == 1  & (type=="vnCDR_wns" | type=="vnCDR+VD_wns") & shots == 0')
# df_plot_1 = density_matrix_test_abs_high.query('nlsp == 1  & (type=="vnCDR_wns" | type=="vnCDR+VD_wns") & shots == 1')
# df_plot_2 = density_matrix_test_abs_high.query('nlsp == 1  & (type=="vnCDR_wns" | type=="vnCDR+VD_wns") & shots == 2')
# df_plot_3 = density_matrix_test_abs_high.query('nlsp == 1  & (type=="vnCDR_wns" | type=="vnCDR+VD_wns") & shots == 3')

# fig, axes = plt.subplots(4,4,figsize=(15,20),sharey=True)
# axes=axes.flatten()

# sns.lineplot(data=df_plot_0.query('Qubits == 4'),x='copies',y='abs_error',hue='type',ax=axes[0],estimator='mean',err_style=None)
# sns.lineplot(data=df_plot_0.query('Qubits == 6'),x='copies',y='abs_error',hue='type',ax=axes[1],estimator='mean',err_style=None)
# sns.lineplot(data=df_plot_0.query('Qubits == 8'),x='copies',y='abs_error',hue='type',ax=axes[2],estimator='mean',err_style=None)
# sns.lineplot(data=df_plot_0.query('Qubits ==10'),x='copies',y='abs_error',hue='type',ax=axes[3],estimator='mean',err_style=None)
# sns.lineplot(data=df_plot_1.query('Qubits == 4'),x='copies',y='abs_error',hue='type',ax=axes[3+1],estimator='mean',err_style=None)
# sns.lineplot(data=df_plot_1.query('Qubits == 6'),x='copies',y='abs_error',hue='type',ax=axes[4+1],estimator='mean',err_style=None)
# sns.lineplot(data=df_plot_1.query('Qubits == 8'),x='copies',y='abs_error',hue='type',ax=axes[5+1],estimator='mean',err_style=None)
# sns.lineplot(data=df_plot_1.query('Qubits ==10'),x='copies',y='abs_error',hue='type',ax=axes[5+2],estimator='mean',err_style=None)
# sns.lineplot(data=df_plot_2.query('Qubits == 4'),x='copies',y='abs_error',hue='type',ax=axes[6+2],estimator='mean',err_style=None)
# sns.lineplot(data=df_plot_2.query('Qubits == 6'),x='copies',y='abs_error',hue='type',ax=axes[7+2],estimator='mean',err_style=None)
# sns.lineplot(data=df_plot_2.query('Qubits == 8'),x='copies',y='abs_error',hue='type',ax=axes[8+2],estimator='mean',err_style=None)
# sns.lineplot(data=df_plot_2.query('Qubits ==10'),x='copies',y='abs_error',hue='type',ax=axes[8+3],estimator='mean',err_style=None)
# sns.lineplot(data=df_plot_3.query('Qubits == 4'),x='copies',y='abs_error',hue='type',ax=axes[9+3],estimator='mean',err_style=None)
# sns.lineplot(data=df_plot_3.query('Qubits ==6'),x='copies',y='abs_error',hue='type',ax=axes[10+3],estimator='mean',err_style=None)
# sns.lineplot(data=df_plot_3.query('Qubits ==8'),x='copies',y='abs_error',hue='type',ax=axes[11+3],estimator='mean',err_style=None)
# sns.lineplot(data=df_plot_3.query('Qubits==10'),x='copies',y='abs_error',hue='type',ax=axes[11+4],estimator='mean',err_style=None)

# for i in axes:
#     i.set_yscale('log')
# axes[0].set_title('4Q - inf shots')
# axes[1].set_title('6Q - inf shots')
# axes[2].set_title('8Q - inf shots')
# axes[3].set_title('10Q - inf shots')
# axes[3+1].set_title('10K shots')
# axes[4+1].set_title('10K  shots')
# axes[5+1].set_title('10K shots')
# axes[5+1+1].set_title('10K shots')
# axes[6+1+1].set_title('100K shots')
# axes[7+1+1].set_title('100K  shots')
# axes[8+1+1].set_title('100K shots')
# axes[8+1+1+1].set_title('100K shots')
# axes[6+3+1+1].set_title('1000K shots')
# axes[7+3+1+1].set_title('1000K  shots')
# axes[8+3+1+1].set_title('1000K shots')
# axes[8+3+1+1+1].set_title('1000K shots')

# # %% 
# fig, axes = plt.subplots(2,3,figsize=(15,10),sharey=True)
# axes=axes.flatten()

# sns.lineplot(data=df_plot_high_s.query('Qubits == 4'),x='copies',y='abs_error',hue='type',ax=axes[1-1],estimator='mean',err_style=None,alpha=0.3)
# sns.lineplot(data=df_plot_high_s.query('Qubits == 6'),x='copies',y='abs_error',hue='type',ax=axes[2-1],estimator='mean',err_style=None,alpha=0.3)
# sns.lineplot(data=df_plot_high_s.query('Qubits == 8'),x='copies',y='abs_error',hue='type',ax=axes[3-1],estimator='mean',err_style=None,alpha=0.3)
# sns.lineplot(data=df_plot_low_s.query('Qubits == 4'),x='copies',y='abs_error',hue='type',ax=axes[3],estimator='mean',err_style=None,alpha=0.3)
# sns.lineplot(data=df_plot_low_s.query('Qubits == 6'),x='copies',y='abs_error',hue='type',ax=axes[4],estimator='mean',err_style=None,alpha=0.3)
# sns.lineplot(data=df_plot_low_s.query('Qubits == 8'),x='copies',y='abs_error',hue='type',ax=axes[5],estimator='mean',err_style=None,alpha=0.3)

# for i in axes:
#     i.set_yscale('log')
# axes[0].set_title('HIGH - 4Q - 100K shots')
# axes[1].set_title('HIGH - 6Q - 100K shots')
# axes[2].set_title('HIGH - 8Q - 100K shots')
# axes[3].set_title('LOW - 4Q - 100K shots')
# axes[4].set_title('LOW - 6Q - 100K shots')
# axes[5].set_title('LOW - 8Q - 100K shots')


# # %%
# df_plot_high = density_matrix_test_abs_high.query('nlsp == 1  & (type=="vnCDR_wns" | type=="vnCDR+VD_wns") & shots == 0 & copies == 6')
# df_plot_low = density_matrix_test_abs_low.query('nlsp == 1  & (type=="vnCDR_wns" | type=="vnCDR+VD_wns") & shots == 0 & copies == 6')
# df_plot_high_s = density_matrix_test_abs_high.query('nlsp == 1  & (type=="vnCDR_wns" | type=="vnCDR+VD_wns") & shots == 1 & copies == 6')
# df_plot_low_s = density_matrix_test_abs_low.query('nlsp == 1  & (type=="vnCDR_wns" | type=="vnCDR+VD_wns") & shots == 1 & copies == 6')

# fig, axes = plt.subplots(2,3,figsize=(15,10))
# axes = axes.flatten()

# sns.boxplot(data=df_plot_high.query('Qubits == 4'),x='type',y='abs_error',ax=axes[1-1])
# sns.boxplot(data=df_plot_high.query('Qubits == 6'),x='type',y='abs_error',ax=axes[2-1])
# sns.boxplot(data=df_plot_high.query('Qubits == 8'),x='type',y='abs_error',ax=axes[3-1])
# sns.boxplot(data=df_plot_low.query('Qubits == 4'),x='type',y='abs_error',ax=axes[3])
# sns.boxplot(data=df_plot_low.query('Qubits == 6'),x='type',y='abs_error',ax=axes[4])
# sns.boxplot(data=df_plot_low.query('Qubits == 8'),x='type',y='abs_error',ax=axes[5])

# for i in axes:
#     i.set_yscale('log')
# axes[0].set_title('HIGH - 4Q - inf shots')
# axes[1].set_title('HIGH - 6Q - inf shots')
# axes[2].set_title('HIGH - 8Q - inf shots')
# axes[3].set_title('LOW - 4Q - inf shots')
# axes[4].set_title('LOW - 6Q - inf shots')
# axes[5].set_title('LOW - 8Q - inf shots')


# # fig, axes = plt.subplots(2,3,figsize=(15,10),sharey=True)
# # axes=axes.flatten()

# # sns.boxplot(data=df_plot_high_s.query('Qubits == 4'),x='type',y='abs_error',ax=axes[1-1])
# # sns.boxplot(data=df_plot_high_s.query('Qubits == 6'),x='type',y='abs_error',ax=axes[2-1])
# # sns.boxplot(data=df_plot_high_s.query('Qubits == 8'),x='type',y='abs_error',ax=axes[3-1])
# # sns.boxplot(data=df_plot_low_s.query('Qubits == 4'),x='type',y='abs_error',ax=axes[3])
# # sns.boxplot(data=df_plot_low_s.query('Qubits == 6'),x='type',y='abs_error',ax=axes[4])
# # sns.boxplot(data=df_plot_low_s.query('Qubits == 8'),x='type',y='abs_error',ax=axes[5])

# # for i in axes:
# #     i.set_yscale('log')
# #     i.set_ylim([10**(-4),10**(-2)])
# # axes[0].set_title('HIGH - 4Q - 100K shots')
# # axes[1].set_title('HIGH - 6Q - 100K shots')
# # axes[2].set_title('HIGH - 8Q - 100K shots')
# # axes[3].set_title('LOW - 4Q - 100K shots')
# # axes[4].set_title('LOW - 6Q - 100K shots')
# # axes[5].set_title('LOW - 8Q - 100K shots')

# # %%

# rescaling_high = density_matrix_test_rescaling_high.query('nlsp == 1  & (type=="VD_wns" ) & shots == 0  ')
# rescaling_high_s = density_matrix_test_rescaling_high.query('nlsp == 1  & (type=="VD_wns" ) & shots == 1 ')
# rescaling_low = density_matrix_test_rescaling_low.query('nlsp == 1  & (type=="VD_wns") & shots == 0 ')
# rescaling_low_s = density_matrix_test_rescaling_low.query('nlsp == 1  & (type=="VD_wns" ) & shots == 1 ')

# fig, axes = plt.subplots(2,2,figsize=(15,10),sharey=True)
# axes=axes.flatten()

# sns.lineplot(data=rescaling_high,x='Qubits',y='error rescaling factor',hue='copies',ax=axes[1-1],estimator='mean',err_style=None)
# sns.lineplot(data=rescaling_high_s,x='Qubits',y='error rescaling factor',hue='copies',ax=axes[2-1],estimator='mean',err_style=None)
# sns.lineplot(data=rescaling_low,x='Qubits',y='error rescaling factor',hue='copies',ax=axes[3-1],estimator='mean',err_style=None)
# sns.lineplot(data=rescaling_low_s,x='Qubits',y='error rescaling factor',hue='copies',ax=axes[3],estimator='mean',err_style=None)

# for i in axes:
#     i.set_yscale('log')
# axes[0].set_title('HIGH  inf shots')
# axes[1].set_title('HIGH - 100K shots')
# axes[2].set_title('LOW - inf shots')
# axes[3].set_title('LOW - 100K shots')

# # %%

# # %%
# simple_abs['Qubits']*simple_abs['depth']
# # %%
# simple_abs['qd']= simple_abs['Qubits']*simple_abs['depth']

# %%

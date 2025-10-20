# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 15:58:46 2024

@author: BioCraze
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import scipy.stats as stats
from statistics import mean, stdev
import numpy as np
from scipy import stats

#Figures for non-viral transfection project

#gcamp_summary = 
#transfected_count = 

fig,axes =plt.subplots(nrows=2, ncols=3, figsize = (25,10), layout="constrained")
palette = ["#95beff","#00bac6"]
sns.set_palette(palette)

ax1 = sns.barplot(x = "transfection_technique", y = "fluo per tectal slice", hue = "Trial", data=gcamp_summary, capsize=0.2, ax=axes[0,0])
ax1.set_ylim(top=max(gcamp_summary["fluo per tectal slice"])*1.2)
sns.stripplot(x = "transfection_technique", y = "fluo per tectal slice", data = gcamp_summary, color="k", size=5, hue="Trial", dodge=True, legend=False, ax=ax1)
ax2=sns.barplot(x = "transfection_technique", y = "fluorescence per ROI", hue = "Trial", data=gcamp_summary, capsize=0.2, ax=axes[0,1])
ax2.set_ylim(top=max(gcamp_summary["fluorescence per ROI"])*1.2)
sns.stripplot(x = "transfection_technique", y = "fluorescence per ROI", data = gcamp_summary, color="k", size=5, hue="Trial", dodge=True, legend=False, ax=ax2)
ax3=sns.barplot(x = "transfection_technique", y = "cell to slice fluo ratio", hue = "Trial", data=gcamp_summary, capsize=0.2, ax=axes[0,2])
ax3.set_ylim(top=max(gcamp_summary["cell to slice fluo ratio"])*1.2)
sns.stripplot(x = "transfection_technique", y = "cell to slice fluo ratio", data = gcamp_summary, color="k", size=5, hue="Trial", dodge=True, legend=False, ax=ax3)
ax4=sns.barplot(x = "transfection_technique", y = "cell count", hue = "Trial", data=gcamp_summary, capsize=0.2, ax=axes[1,0])
ax4.set_ylim(top=max(gcamp_summary["cell count"])*1.2)
sns.stripplot(x = "transfection_technique", y = "cell count", data = gcamp_summary, color="k", size=5, hue="Trial", dodge=True, legend=False, ax=ax4)
ax5=sns.barplot(x = "transfection_technique", y = "per cell fluo", hue = "Trial", data=cell_fluo, capsize=0.2, ax=axes[1,1])
#ax5.set_ylim(top=max(cell_fluo["per cell fluo"])*1.2)
#sns.stripplot(x = "transfection_technique", y = "per cell fluo", data = cell_fluo, color="k", size=3, hue="Trial", dodge =True,ax=ax5)
ax6=sns.barplot(x = "Transfection technique", y = "transfected cell count", data=transfected_count, capsize=0.2, palette=["#00a258","#827000","#ad0f00"], ax=axes[1,2])
ax6.set_ylim(top=max(transfected_count["transfected cell count"])*1.2)
sns.stripplot(x = "Transfection technique", y = "transfected cell count", data = transfected_count, color="k", size=5, legend=False, ax=ax6)
plt.show()
#------------------------------------------------------------------------------------------
colour='white'
mpl.rcParams["text.color"]=colour
mpl.rcParams["axes.labelcolor"]=colour
mpl.rcParams["xtick.color"]=colour
mpl.rcParams["ytick.color"]=colour
mpl.rcParams["axes.facecolor"]='k'
mpl.rcParams["axes.edgecolor"]=colour
mpl.rcParams["figure.edgecolor"]=colour
mpl.rcParams["figure.facecolor"]='k'
mpl.rcParams["font.size"]=15
palette = ["#1be6ec", "#e51bec"]

#Figures for glia plasticity
#GLIA DATA
    #no training
fig,axes =plt.subplots(nrows=2, ncols=4, figsize = (25,10), layout="constrained")
sns.set_palette(palette)

responsive_glia_data= pd.read_csv("E:/glia projects/plasticity/summaries/glia_training_glia_responding_cell_numbers.csv").sort_values(by=['time group']).replace(to_replace='min99',value='min100')
responsive_glia_data=responsive_glia_data.loc[~((responsive_glia_data['time group']=='min30') | (pd.isna(responsive_glia_data['time group'])))]
responsive_glia_ctrl = sns.pointplot(x = "time group", y = "Corrected cell count", hue = "treatment", data=responsive_glia_data[responsive_glia_data['treatment'].str.contains('without training')], errorbar='se', capsize=0.2, ax=axes[0,0])
responsive_glia_ctrl.set_ylim(top=max(responsive_glia_data["Corrected cell count"])*1.2, bottom=min(responsive_glia_data["Corrected cell count"])*0.6)
responsive_glia_ctrl.set_title("Active glia count - no training")
sns.stripplot(x = "time group", y = "Corrected cell count", data = responsive_glia_data[responsive_glia_data['treatment'].str.contains('without training')], size=5, hue="treatment", dodge=False, legend=False, ax=responsive_glia_ctrl)

responsive_glia_exp= sns.pointplot(x = "time group", y = "Corrected cell count", hue = "treatment", data=responsive_glia_data[responsive_glia_data['treatment'].str.contains('with training')], errorbar='se', capsize=0.2, ax=axes[1,0])
responsive_glia_exp.set_ylim(top=max(responsive_glia_data["Corrected cell count"])*1.2, bottom=min(responsive_glia_data["Corrected cell count"])*0.6)
responsive_glia_exp.set_title("Active glia count - training")
sns.stripplot(x = "time group", y = "Corrected cell count", data = responsive_glia_data[responsive_glia_data['treatment'].str.contains('with training')], size=5, hue="treatment", dodge=False, legend=False, ax=responsive_glia_exp)


glia_transients_count_data= pd.read_csv("E:/glia projects/plasticity/summaries/glia_training_glia_response_numbers.csv").sort_values(by=['time group']).replace(to_replace='min99',value='min100')
glia_transients_count_data=glia_transients_count_data.loc[~((glia_transients_count_data['time group']=='min30') | (pd.isna(glia_transients_count_data['time group'])))]
glia_transients_count_ctrl = sns.pointplot(x = "time group", y = "normalized response count", hue = "treatment", data=glia_transients_count_data[glia_transients_count_data['treatment'].str.contains('without training')], errorbar='se', capsize=0.2, ax=axes[0,1])
glia_transients_count_ctrl.set_ylim(top=max(glia_transients_count_data["normalized response count"])*1.2, bottom=min(glia_transients_count_data["normalized response count"])*0.6)
glia_transients_count_ctrl.set_title("Glia transients counts - no training")
sns.stripplot(x = "time group", y = "normalized response count", data = glia_transients_count_data[glia_transients_count_data['treatment'].str.contains('without training')], size=5, hue="treatment", dodge=False, legend=False, ax=glia_transients_count_ctrl)
	
	#With training
glia_transients_count_exp= sns.pointplot(x = "time group", y = "normalized response count", hue = "treatment", data=glia_transients_count_data[glia_transients_count_data['treatment'].str.contains('with training')], errorbar='se', capsize=0.2, ax=axes[1,1])
glia_transients_count_exp.set_ylim(top=max(glia_transients_count_data["normalized response count"])*1.2, bottom=min(glia_transients_count_data["normalized response count"])*0.6)
glia_transients_count_exp.set_title("Glia transients counts - training")
sns.stripplot(x = "time group", y = "normalized response count", data = glia_transients_count_data[glia_transients_count_data['treatment'].str.contains('with training')], size=5, hue="treatment", dodge=False, legend=False, ax=glia_transients_count_exp)

glia_transients_data= pd.read_csv("E:/glia projects/plasticity/summaries/glia_training_glia_summary_by_cell.csv").sort_values(by=['time group']).replace(to_replace='min99',value='min100')
glia_transients_data=glia_transients_data.loc[~((glia_transients_data['time group']=='min30') | (pd.isna(glia_transients_data['time group'])))]
#glia_transients_amp_ctrl = sns.pointplot(x = "time group", y = "peaks", hue = "treatment", data=glia_transients_data[glia_transients_data['treatment'].str.contains('without training')], errorbar='se', capsize=0.2, ax=axes[0,2])
#glia_transients_amp_ctrl.set_ylim(top=max(glia_transients_data["peaks"])*1.2, bottom=min(glia_transients_data["peaks"])*0.6)
#glia_transients_amp_ctrl.set_title("Glia peak amplitude - no training")
#sns.stripplot(x = "time group", y = "peaks", data = glia_transients_data[glia_transients_data['treatment'].str.contains('without training')], size=5, hue="treatment", dodge=False, legend=False, ax=glia_transients_amp_ctrl)

glia_transients_amp_exp= sns.pointplot(x = "time group", y = "peaks", hue = "treatment", data=glia_transients_data[glia_transients_data['treatment'].str.contains('with training')], errorbar='se', capsize=0.2, ax=axes[1,2])
glia_transients_amp_exp.set_ylim(top=max(glia_transients_data["peaks"])*1.2, bottom=min(glia_transients_data["peaks"])-min(glia_transients_data["peaks"])*1.2)
glia_transients_amp_exp.set_title("Glia peak amplitude - training")
sns.stripplot(x = "time group", y = "peaks", data = glia_transients_data[glia_transients_data['treatment'].str.contains('with training')], size=5, hue="treatment", dodge=False, legend=False, ax=glia_transients_amp_exp)

glia_transients_data= pd.read_csv("E:/glia projects/plasticity/summaries/glia_training_glia_summary_by_cell.csv").sort_values(by=['time group']).replace(to_replace='min99',value='min100')
glia_transients_data=glia_transients_data.loc[~((glia_transients_data['time group']=='min30') | (pd.isna(glia_transients_data['time group'])))]
#glia_transients_amp_ctrl = sns.pointplot(x = "time group", y = "area", hue = "treatment", data=glia_transients_data[glia_transients_data['treatment'].str.contains('without training')], errorbar='se', capsize=0.2, ax=axes[0,3])
#glia_transients_amp_ctrl.set_ylim(top=max(glia_transients_data["area"])*1.2, bottom=min(glia_transients_data["area"])*0.6)
#glia_transients_amp_ctrl.set_title("Glia auc - no training")
#sns.stripplot(x = "time group", y = "area", data = glia_transients_data[glia_transients_data['treatment'].str.contains('without training')], size=5, hue="treatment", dodge=False, legend=False, ax=glia_transients_amp_ctrl)

glia_transients_amp_exp= sns.pointplot(x = "time group", y = "area", hue = "treatment", data=glia_transients_data[glia_transients_data['treatment'].str.contains('with training')], errorbar='se', capsize=0.2, ax=axes[1,3])
glia_transients_amp_exp.set_ylim(top=max(glia_transients_data["area"])*1.2, bottom=min(glia_transients_data["area"])-min(glia_transients_data["area"])*1.2)
glia_transients_amp_exp.set_title("Glia auc - training")
sns.stripplot(x = "time group", y = "area", data = glia_transients_data[glia_transients_data['treatment'].str.contains('with training')], size=5, hue="treatment", dodge=False, legend=False, ax=glia_transients_amp_exp)

plt.show()


#NEURON DATA
neuron_response_data= pd.read_csv("E:/glia projects/plasticity/summaries/glia_training_old_neuron_summary.csv").sort_values(by=['time group'])
neuron_response_data=neuron_response_data.loc[~((neuron_response_data['time group']==30) | (pd.isna(neuron_response_data['time group'])))]
neuron_response_data_by_anml = neuron_response_data[["Animal","corrected peaks","corrected area","treatment","time group"]].groupby(['Animal','treatment','time group']).mean().reset_index()
neuron_response_data_by_cell = neuron_response_data[["Animal","corrected peaks","corrected area","treatment","time group", 'trace']].groupby(['Animal','trace','treatment','time group']).mean().reset_index()

    #No training
data_notrain_by_anml=neuron_response_data_by_anml[neuron_response_data_by_anml['treatment'].str.contains('without training')]
data_notrain_by_cell=neuron_response_data_by_cell[neuron_response_data_by_cell['treatment'].str.contains('without training')]

fig,axes =plt.subplots(nrows=1, ncols=1, figsize = (20,10), layout="constrained",)
sns.set_palette(palette)

neuron_response_amp_ctrl = sns.pointplot(data=data_notrain_by_cell, x = "time group", y = "corrected peaks", hue = "treatment", errorbar='se', capsize=0.1, err_kws={'color':'white','linewidth':2}, ax=axes)
neuron_response_amp_ctrl.set_ylim(top=1.5, bottom=0.5)
neuron_response_amp_ctrl.set_title("Neuron peak amplitude - no training")
neuron_response_amp_ctrl.set_ylabel("Normalized response amplitude (F/F0)")
neuron_response_amp_ctrl.set_xlabel("Time (min)")
sns.stripplot(data = data_notrain_by_cell, x = "time group", y = "corrected peaks", size=5, hue="treatment", dodge=False, legend=False, ax=neuron_response_amp_ctrl, alpha=0.4)


fig,axes =plt.subplots(nrows=1, ncols=1, figsize = (20,10), layout="constrained")
neuron_response_auc_ctrl = sns.pointplot(x = "time group", y = "corrected area", hue = "treatment", data=data_notrain_by_cell, errorbar='se', capsize=0.1, err_kws={'color':'white','linewidth':2}, ax=axes)
neuron_response_auc_ctrl.set_ylim(top=1.5, bottom=0.5)
neuron_response_auc_ctrl.set_title("Neuron auc - no training")
neuron_response_auc_ctrl.set_ylabel("Normalized response amplitude (F/F0)")
neuron_response_auc_ctrl.set_xlabel("Time (min)")
sns.stripplot(x = "time group", y = "corrected area", data = data_notrain_by_cell, size=5, hue="treatment", dodge=False, legend=False, ax=neuron_response_auc_ctrl, alpha=0.4)

    
    #With training
data_training_by_anml = neuron_response_data_by_anml[neuron_response_data_by_anml['treatment'].str.contains('with training')]
data_training_by_cell = neuron_response_data_by_cell[neuron_response_data_by_cell['treatment'].str.contains('with training')]
data_training = neuron_response_data[neuron_response_data['treatment'].str.contains('with training')]
fig,axes =plt.subplots(nrows=1, ncols=1, figsize = (20,10), layout="constrained")
sns.set_palette(palette)

neuron_response_amp_exp= sns.pointplot(x = "time group", y = "corrected peaks", hue = "treatment", data=data_training, errorbar='se', capsize=0.1, err_kws={'color':'white','linewidth':2}, ax=axes)
neuron_response_amp_exp.set_ylim(top=1.2, bottom=0.6)
neuron_response_amp_exp.set_title("Neuron peak amplitude - training")
neuron_response_amp_exp.set_ylabel("Normalized response amplitude (F/F0)")
neuron_response_amp_exp.set_xlabel("Time (min)")
sns.stripplot(x = "time group", y = "corrected peaks", data = data_training, size=5, hue="treatment", dodge=False, legend=False, ax=neuron_response_amp_exp, alpha=0.4)

fig,axes =plt.subplots(nrows=1, ncols=1, figsize = (20,10), layout="constrained")
neuron_response_auc_exp= sns.pointplot(x = "time group", y = "corrected area", hue = "treatment", data=data_training, errorbar='se', capsize=0.1, err_kws={'color':'white','linewidth':2}, ax=axes)
neuron_response_auc_exp.set_ylim(top=1.2, bottom=0.6)
neuron_response_auc_exp.set_title("Neuron auc - training")
neuron_response_auc_exp.set_ylabel("Normalized response amplitude (F/F0)")
neuron_response_auc_exp.set_xlabel("Time (min)")
sns.stripplot(x = "time group", y = "corrected area", data = data_training, size=5, hue="treatment", dodge=False, legend=False, ax=neuron_response_auc_exp, alpha=0.4)


#NEUROPIL DATA
neuropil_data_norm = pd.read_csv("E:/glia projects/plasticity/summaries/glia training neuropil response summary.csv").sort_values(by=['time group'])
neuropil_data_norm=neuropil_data_norm.loc[~((neuropil_data_norm['time group']==30) | (pd.isna(neuropil_data_norm['time group'])))]
neuropil_data_norm_by_anml = neuropil_data_norm[["animal","corrected peaks","corrected area","treatment","time group"]].groupby(['animal','treatment','time group']).mean().reset_index()
neuropil_data_raw = pd.read_csv("E:/glia projects/plasticity/summaries/neuropil_au_peaks.csv").sort_values(by=['time group'])
neuropil_data_raw=neuropil_data_raw.loc[~((neuropil_data_raw['time group']==30) | (pd.isna(neuropil_data_raw['time group'])))]
neuropil_data_raw_by_anml = neuropil_data_raw[["animal","Normalized au peaks","treatment","time group"]].groupby(['animal','treatment','time group']).mean().reset_index()

    #without training
        #delta f
neuropil_norm_anml_ctrl = neuropil_data_norm_by_anml[neuropil_data_norm_by_anml['treatment'].str.contains('without training')]
fig,axes =plt.subplots(nrows=1, ncols=1, figsize = (20,10), layout="constrained")
neuropil_norm_auc_ctrl= sns.pointplot(x = "time group", y = "corrected area", hue = "treatment", data=neuropil_norm_anml_ctrl, errorbar='se', capsize=0.1, err_kws={'color':'white','linewidth':2}, ax=axes)
#neuron_response_auc_exp.set_ylim(top=1.5, bottom=0.5)
neuropil_norm_auc_ctrl.set_title("Neuropil auc - without training")
neuropil_norm_auc_ctrl.set_ylabel("Normalized response AUC (F/F0)")
neuropil_norm_auc_ctrl.set_xlabel("Time (min)")
sns.stripplot(x = "time group", y = "corrected area", data = neuropil_norm_anml_ctrl, size=5, hue="treatment", dodge=False, legend=False, ax=neuropil_norm_auc_ctrl, alpha=1)
        
        #raw
neuropil_raw_anml_ctrl = neuropil_data_raw_by_anml[neuropil_data_raw_by_anml['treatment'].str.contains('without training')]
fig,axes =plt.subplots(nrows=1, ncols=1, figsize = (20,10), layout="constrained")
neuropil_raw_amp_ctrl= sns.pointplot(x = "time group", y = "Normalized au peaks", hue = "treatment", data=neuropil_raw_anml_ctrl, errorbar='se', capsize=0.1, err_kws={'color':'white','linewidth':2}, ax=axes)
#neuron_response_auc_exp.set_ylim(top=1.5, bottom=0.5)
neuropil_raw_amp_ctrl.set_title("Neuropil amplitude - without training")
neuropil_raw_amp_ctrl.set_ylabel("Normalized response amplitude (a.u.)")
neuropil_raw_amp_ctrl.set_xlabel("Time (min)")
sns.stripplot(x = "time group", y = "Normalized au peaks", data = neuropil_raw_anml_ctrl, size=5, hue="treatment", dodge=False, legend=False, ax=neuropil_raw_amp_ctrl, alpha=1)


    #with training
        #delta f
neuropil_norm_anml_training = neuropil_data_norm_by_anml[neuropil_data_norm_by_anml['treatment'].str.contains('with training')]
fig,axes =plt.subplots(nrows=1, ncols=1, figsize = (20,10), layout="constrained")
neuropil_response_auc_exp= sns.pointplot(x = "time group", y = "corrected area", hue = "treatment", data=neuropil_norm_anml_training, errorbar='se', capsize=0.1, ax=axes)
#neuron_response_auc_exp.set_ylim(top=1.5, bottom=0.5)
neuropil_response_auc_exp.set_title("Neuropil auc - with training")
neuropil_response_auc_exp.set_ylabel("Normalized response AUC (F/F0)")
neuropil_response_auc_exp.set_xlabel("Time (min)")
sns.stripplot(x = "time group", y = "corrected area", data = neuropil_norm_anml_training, size=5, hue="treatment", dodge=False, legend=False, ax=neuropil_response_auc_exp, alpha=1)
        
        #raw
neuropil_raw_anml_exp = neuropil_data_raw_by_anml[neuropil_data_raw_by_anml['treatment'].str.contains('with training')]
fig,axes =plt.subplots(nrows=1, ncols=1, figsize = (20,10), layout="constrained")
neuropil_raw_amp_exp= sns.pointplot(x = "time group", y = "Normalized au peaks", hue = "treatment", data=neuropil_raw_anml_exp, errorbar='se', capsize=0.1, ax=axes)
#neuron_response_auc_exp.set_ylim(top=1.5, bottom=0.5)
neuropil_raw_amp_exp.set_title("Neuropil amplitude - with training")
neuropil_raw_amp_exp.set_ylabel("Normalized response amplitude (a.u.)")
neuropil_raw_amp_exp.set_xlabel("Time (min)")
sns.stripplot(x = "time group", y = "Normalized au peaks", data = neuropil_raw_anml_exp, size=5, hue="treatment", dodge=False, legend=False, ax=neuropil_raw_amp_exp, alpha=1)


#CELL LOCATION MATCHING ------------------------------------------------------------------------------------------
fig,ax =plt.subplots(figsize = (25,25), layout="constrained")
ax = sns.heatmap(fltrd_dst_matrix, vmin=np.min(fltrd_dst_matrix), vmax=np.min(fltrd_dst_matrix)*100 ,center=np.min(fltrd_dst_matrix)*10, cbar=False)
ax.set_ylabel("cell id in current time point")
ax.set_xlabel("cell id in next time point")
ax.set_title("Matching cells")

#------------------------------------------------------------------------------------------------------------------------------
#FIGURES FOR LOCAL FIELD POTENTIALS

field_data = pd.read_csv('E:/glia projects/field recordings/field_potential_data.csv')
field_data=field_data[(field_data["series time"]>=0)&(field_data["series time"]<=47.5)].sort_values(by=["treatment"],ascending=False)
field_stats = {}
times=np.arange(0,50,2.5)

#stats
for a in pd.unique(field_data['treatment']):
    for b in pd.unique(field_data['treatment']):
        if a!=b:

            field_stats[a+"vs"+b]={"pvalue":[],"test":[]}
            for time in times:
                sa=field_data['normalized peak amp'][(field_data['series time']==time) & (field_data['treatment']==a)]
                sb=field_data['normalized peak amp'][(field_data['series time']==time) & (field_data['treatment']==b)]
                if stats.shapiro(sa,nan_policy='omit').pvalue>0.01 and stats.shapiro(sb,nan_policy='omit').pvalue>0.01:
                    if stats.levene(sa,sb,axis=0,nan_policy='omit').pvalue>0.01:
                        field_stats[a+"vs"+b]['pvalue'].append(stats.ttest_ind(sa,sb,equal_var=True,nan_policy='omit').pvalue)
                        field_stats[a+"vs"+b]['test'].append("ttest_equal_var")
                    else:
                        field_stats[a+"vs"+b]['pvalue'].append(stats.ttest_ind(sa,sb,equal_var=False,nan_policy='omit').pvalue)
                        field_stats[a+"vs"+b]['test'].append("ttest_unequal_var")
                else:
                    field_stats[a+"vs"+b]['pvalue'].append(stats.mannwhitneyu(sa, sb, nan_policy='omit').pvalue)
                    field_stats[a+"vs"+b]['test'].append("mwu")


field_stats_df = pd.DataFrame(field_stats).T.iloc[[0,1,3]]

corrected_pvalues = stats.false_discovery_control(list(field_stats_df["pvalue"]),axis=None)
corrected_pvalues_list=[list(corrected_pvalues[0:20]),list(corrected_pvalues[20:40]),list(corrected_pvalues[40:60])]
field_stats_df["corrected pvalues"]=corrected_pvalues_list
reconfig={}
for treatment in field_stats_df.index:
    for i in field_stats_df.columns:
        reconfig[treatment,i]=field_stats_df.loc[treatment,i]
field_stats_df_reconfig=pd.DataFrame(reconfig).T
field_stats_df_reconfig.to_csv("E:/glia projects/field recordings/field_potential_stats.csv")

#peak amplitude
fig,axes =plt.subplots(nrows=1, ncols=1, figsize = (20,10), layout="constrained",)
sns.set_palette(palette)

field_peak_amp = sns.pointplot(data=field_data, x = "series time", y = "normalized peak amp", hue = "treatment", errorbar='se', capsize=0.1, err_kws={'color':'white','linewidth':2}, ax=axes)
#field_peak_amp.set_ylim(top=2, bottom=0)
field_peak_amp.set_title("LFP normalized peak amplitude")
field_peak_amp.set_ylabel("Normalized amplitude (F/F0)")
field_peak_amp.set_xlabel("Time (min)")
sns.stripplot(data = field_data, x = "series time", y = "normalized peak amp", size=5, hue="treatment", dodge=False, legend=False, ax=field_peak_amp, alpha=0.4)

#absolute amplitude
fig,axes =plt.subplots(nrows=1, ncols=1, figsize = (20,10), layout="constrained",)
sns.set_palette(palette)

field_abs_amp = sns.pointplot(data=field_data, x = "series time", y = "normalized absolute amp", hue = "treatment", errorbar='se', capsize=0.1, err_kws={'color':'white','linewidth':2}, ax=axes)
#field_peak_amp.set_ylim(top=2, bottom=0)
field_abs_amp.set_title("LFP normalized absolute amplitude")
field_abs_amp.set_ylabel("Normalized amplitude (F/F0)")
field_abs_amp.set_xlabel("Time (min)")
sns.stripplot(data = field_data, x = "series time", y = "normalized absolute amp", size=5, hue="treatment", dodge=False, legend=False, ax=field_abs_amp, alpha=0.4)

#before and after barchart for absolute amp
pre_lfp_abs_amp = field_data[(field_data["series time"]==0)|(field_data["series time"]==2.5)][["normalized absolute amp","treatment"]]
pre_lfp_abs_amp["time"]="t0"
post_lfp_abs_amp = field_data[(field_data["series time"]==45.0)|(field_data["series time"]==47.5)][["normalized absolute amp","treatment"]]
post_lfp_abs_amp["time"]="t50"

pre_post_lfp_abs_amp = pd.concat([pre_lfp_abs_amp,post_lfp_abs_amp]).sort_values(by=['treatment'],ascending=False)

fig,axes =plt.subplots(nrows=1, ncols=1, figsize = (10,10), layout="constrained",)
sns.set_palette(palette)

field_abs_amp = sns.barplot(data=pre_post_lfp_abs_amp, x = "time", y = "normalized absolute amp", hue = "treatment", errorbar='se', capsize=0.1, err_kws={'color':'white','linewidth':2}, ax=axes)
#field_peak_amp.set_ylim(top=2, bottom=0)
field_abs_amp.set_title("LFP normalized absolute amplitude")
field_abs_amp.set_ylabel("Normalized amplitude (F/F0)")
field_abs_amp.set_xlabel("Time (min)")
sns.stripplot(data = pre_post_lfp_abs_amp, x = "time", y = "normalized absolute amp", size=5, hue="treatment", dodge=True, linewidth=1, edgecolor="white",legend=False, ax=field_abs_amp)


#fiber volley
fig,axes =plt.subplots(nrows=1, ncols=1, figsize = (20,10), layout="constrained",)
sns.set_palette(palette)

fv_amp = sns.pointplot(data=field_data, x = "series time", y = "normalized fv amp", hue = "treatment", errorbar='se', capsize=0.1, err_kws={'color':'white','linewidth':2}, ax=axes)
#fv_amp.set_ylim(top=2, bottom=0)
fv_amp.set_title("LFP normalized fiber volley amplitude")
fv_amp.set_ylabel("Normalized amplitude (F/F0)")
fv_amp.set_xlabel("Time (min)")
sns.stripplot(data = field_data, x = "series time", y = "normalized fv amp", size=5, hue="treatment", dodge=False, legend=False, ax=fv_amp, alpha=0.4)







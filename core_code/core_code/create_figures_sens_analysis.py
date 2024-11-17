
# importing packages
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import pandas as pd
import numpy as np

data_out_folder = "../../../../results/Data_out/"
#Load sens_analysis file
## (I have seperated the modelrun_id from the rest of the config file so that I can put it in the .gitignore)
dfile = pd.read_csv(data_out_folder + "peak_jobs_sum_diff_from_ref_fig.csv", header = 0, dtype={'Group':str,'ModelRun':str,'Description':str,'Value':str,'PeakJobs':float,'TextMove':float}, index_col=False)

df_IO = dfile.loc[dfile['Group'] == 'IO']
df_IO = df_IO.reset_index(drop=True)
df_LC = dfile.loc[dfile['Group'] == 'Line Cost']
df_LC = df_LC.reset_index(drop=True)
df_TD = dfile.loc[dfile['Group'] == 'T&D opex']
df_TD = df_TD.reset_index(drop=True)
df_OI = dfile.loc[dfile['Group'] == 'Occ Ind']
df_OI = df_OI.reset_index(drop=True)
df_Sm = dfile.loc[dfile['Group'] == 'Smooth']
df_Sm = df_Sm.reset_index(drop=True)
df_CC = dfile.loc[dfile['Group'] == 'Cost curves']
df_CC = df_CC.reset_index(drop=True)
df_CV = dfile.loc[dfile['Group'] == 'Capex vec']
df_CV = df_CV.reset_index(drop=True)
df_OV = dfile.loc[dfile['Group'] == 'Opex vec']
df_OV = df_OV.reset_index(drop=True)
df_All = dfile.loc[dfile['Group'] == 'All']
df_All = df_All.reset_index(drop=True)

for jobtype in ['PeakJobs', 'SteadyJobs']:

    afig, (ax,ax2) = plt.subplots(2,1,figsize = (10,7),gridspec_kw={'height_ratios': [10, 1]})
    ax.get_yaxis().set_major_formatter(tick.FuncFormatter(lambda x, p: format(int(x), ',')))
    #ymin = df_All.loc[:,jobtype].min()
    #ymax = df_All.loc[:,jobtype].max()
    if jobtype == 'PeakJobs':
        ymin = 400000
        ymax = 900000
    elif jobtype == 'SteadyJobs':
        ymin = -50000
        ymax = 250000
    #ax.set_ylim(ymin*1.1 if ymin < 0 else ymin*0.9,ymax*1.1)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0.5,7.5)

    def plot_sens_range(df, col_num,ax):
        # making a simple plot
        the_min = df.loc[:,jobtype].min()
        the_max = df.loc[:,jobtype].max()
    
        ax.vlines(x = col_num, ymin = the_min, ymax = the_max, colors = 'black')
        size_of_dots = 5 #size of dots for graphic
        ax.scatter([col_num]*len(df.loc[:,jobtype]),df.loc[:,jobtype],color='k',s=size_of_dots)
        for i in range(df.shape[0]):
            #Add text if it is not NAN
            if pd.notna(df.loc[i,'Value']):
                if jobtype == 'PeakJobs':
                    ax.text(x=col_num+0.1,y=float(df.loc[i,jobtype])+float(df.loc[i,'TextMove']),s=str(df.loc[i,'Value']),fontdict=dict(color='black',size=7))
                elif jobtype == 'SteadyJobs':
                    ax.text(x=col_num+0.1,y=float(df.loc[i,jobtype])+float(df.loc[i,'TextMoveSteadyState']),s=str(df.loc[i,'Value']),fontdict=dict(color='black',size=7))

        #plt.text([col_num]*len(df.loc[:,'PeakJobs']),df.loc[:,'PeakJobs'],df.loc[:,'Value'])
    

    plot_sens_range(df_IO, 1,ax)
    plot_sens_range(df_LC, 2,ax)
    plot_sens_range(df_TD, 3,ax)
    plot_sens_range(df_OI, 4,ax)
    #plot_sens_range(df_Sm, 5,ax)
    plot_sens_range(df_CC, 5,ax)
    plot_sens_range(df_CV, 6,ax)
    plot_sens_range(df_OV, 7,ax)
    #plot_sens_range(df_All, 8,ax)

    bars = ['I/O table base year', 'Transmission & Distribution line costs', 'Transmission & Distribution opex costs', 'Occupation-industry employment', 'Technology cost curves', 'Capex cost vectors (30 model runs with random noise)', 'Opex cost vectors (30 model runs with random noise)']#,'All']
    ax.set_xticklabels(ax.get_xticks(),rotation = 90,wrap=True,horizontalalignment='right')
    ax.set_xticks(np.arange(1,8), bars)
    ax.set_xlabel("Sensitivity Analysis")
    if jobtype == 'PeakJobs':
        ax.set_ylabel("Estimated new jobs in 2034")
    elif jobtype == 'SteadyJobs':
        ax.set_ylabel("Estimated new jobs in 2045")

    # make room for x-axis
    afig.subplots_adjust(bottom=0.08)
    ax.grid(True, which='major', axis='y', lw=0.3, alpha=0.6)
    # put a horizontal line on the med value
    if jobtype == 'PeakJobs':
        ax.hlines(y=582345.74, xmin = 0, xmax = 9, colors = 'grey',linestyles='dashed', alpha=0.4)
    else: # SteadyJobs
        ax.hlines(y=97914.66, xmin = 0, xmax = 9, colors = 'grey',linestyles='dashed', alpha=0.4)

    ax2.xaxis.set_visible(False) 
    ax2.yaxis.set_visible(False) 
    ax2.axis('off')

    #save figure
    plt.savefig(data_out_folder + jobtype + '_sens_analysis_fig.png')
    plt.savefig(data_out_folder + jobtype + '_sens_analysis_fig.pdf')
    plt.savefig(data_out_folder + jobtype + '_sens_analysis_fig.jpeg')
    plt.show()


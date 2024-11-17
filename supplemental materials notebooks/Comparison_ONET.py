import pandas as pd
import numpy as np
from matplotlib import pylab as plt
import sklearn as sk

from sklearn.metrics import precision_score, recall_score

path_data_omn = "../data/omn/"
path_data_scenarios = "../results/data_out/results/"
path_data_out = "../results/data_out/corr_green_brown/"
path_fig = "../results/figs/"
file_arch = "occ_archetypes_thresholds_relbase_2034_2038.csv"
#file_arch = "occ_archetypes_thresholds.csv"
file_shock_2044 = "diff_baseline_employment_effects_2044-2045__3_oct_2022.csv"
shock_name = "diff_baseline_employment_effects_"
shock_verison = "__3_oct_2022.csv"

start_year = 2021
end_year = 2044
mark_year = 2035 # year marking the before after is inclusive

####################
# import archetypes
###################
df_shock_and_archetypes = pd.read_csv(path_data_omn + file_arch)

archaetypes = ['Phase_out', 'Permanent_boost', 'Temporary_boost', 'Late_boost']

archaetypes = ['Phase_out_r0.01', 'Permanent_boost_r0.01',\
    'Temporary_boost_r0.01', 'Late_boost_r0.01']

onet_measures = ['Green Enhanced Skills', \
    'Green New & Emerging','Green Increased Demand', 'brown']

df_archetypes = pd.read_csv(path_data_omn + file_arch)
df_archetypes.drop(['total_shock_positive', 'total_shock_negative', 
        'total_shock_opex','total_shock_opex_positive', 
        'total_shock_opex_negative','total_shock_capex',
        'total_shock_capex_positive','total_shock_capex_negative',], axis=1)
# shock variables

#########
# Get O*NET green skills
#########
# Get O*NET green skills
df_onet = pd.read_csv(path_data_omn + "Green_Occupations.csv")
# get until six digits
# NOTE this means we consider an occupation green if one sub classification
# of the occupation is green.
df_onet['O*NET-SOC Code'] = df_onet['O*NET-SOC Code'].str.slice(start=0, stop=7)
# also put those ending with zero for when shocks are not as detailed
df_onet['O*NET-SOC Code broad'] = df_onet['O*NET-SOC Code'].str.slice(start=0, stop=6)
df_onet['O*NET-SOC Code broad'] = df_onet['O*NET-SOC Code broad'] + "0"
# Make dummies of green
df_onet = pd.get_dummies(df_onet,columns=['Green Occupational Category'])

# Make dictionaries for soc codes six and five digit of greeness
dict_code_green_enhanced = dict(zip(df_onet["O*NET-SOC Code"], \
                        df_onet["Green Occupational Category_Green Enhanced Skills"]))
dict_code_green_increase = dict(zip(df_onet["O*NET-SOC Code"], \
                        df_onet["Green Occupational Category_Green Increased Demand"]))
dict_code_green_emerging = dict(zip(df_onet["O*NET-SOC Code"], \
                        df_onet["Green Occupational Category_Green New & Emerging"]))
# add them at a broader level
dict_code_green_enhanced_broad = dict(zip(df_onet["O*NET-SOC Code broad"], \
                        df_onet["Green Occupational Category_Green Enhanced Skills"]))
dict_code_green_increase_broad = dict(zip(df_onet["O*NET-SOC Code broad"], \
                        df_onet["Green Occupational Category_Green Increased Demand"]))
dict_code_green_emerging_broad = dict(zip(df_onet["O*NET-SOC Code broad"], \
                        df_onet["Green Occupational Category_Green New & Emerging"]))
# Merge both of them
dict_code_green_enhanced = {**dict_code_green_enhanced, **dict_code_green_enhanced_broad}
dict_code_green_increase = {**dict_code_green_increase, **dict_code_green_increase_broad}
dict_code_green_emerging = {**dict_code_green_emerging, **dict_code_green_emerging_broad}


#####
# Import Vona Brown
#####
df_brown = pd.read_csv(path_data_omn +  "brown_occ_vona.csv", delimiter=";")
# also put those codes ending with zero for when shocks are not as detailed
df_brown['SOC_2010 broad'] = df_brown['SOC_2010'].str.slice(start=0, \
    stop=6)
df_brown['SOC_2010 broad'] = df_brown['SOC_2010 broad'] + "0"

list_brown = df_brown["SOC_2010"].to_list() + df_brown["SOC_2010 broad"].to_list()


####
# Compare with archetypes
####

df_archetypes['Green Enhanced Skills'] = df_archetypes['O*NET-SOC Code']\
    .map(dict_code_green_enhanced )
df_archetypes['Green Increased Demand'] = df_archetypes['O*NET-SOC Code']\
    .map(dict_code_green_increase )
df_archetypes['Green New & Emerging'] = df_archetypes['O*NET-SOC Code']\
    .map(dict_code_green_emerging )

df_archetypes = df_archetypes.fillna(0)


df_archetypes['brown'] = np.where(df_archetypes["O*NET-SOC Code"]\
    .isin(list_brown), 1, 0)


df_selarch_green = df_archetypes[[ 'Phase_out_r0.01', 'Permanent_boost_r0.01',\
    'Temporary_boost_r0.01', 'Late_boost_r0.01', 'Green Enhanced Skills', \
    'Green New & Emerging','Green Increased Demand', 'brown']]

from scipy.stats import pearsonr
import numpy as np
rho = df_selarch_green.corr()
pval = df_selarch_green.corr(method=lambda x, y: pearsonr(x, y)[1]) - np.eye(*rho.shape)
p = pval.applymap(lambda x: ''.join(['*' for t in [.05, .01, .001] if x<=t]))
rho.round(2).astype(str) + p
print("pearson correlation")
print(rho.round(2).astype(str) + p)

# save rho to csv
rho.to_csv(path_data_out + "rho.csv")
# include stars (from p value) in csv
(rho.round(2).astype(str) + p).to_csv(path_data_out + "rho_stars.csv")
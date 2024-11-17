import pandas as pd
import numpy as np
from matplotlib import pylab as plt
import seaborn as sns
import string
import scipy.stats
import matplotlib.colors
from sklearn.metrics import precision_score, recall_score
#paths
path_data_maria = "../data/omn/"

path_fig = "../Data_out/fig/"
# shocks
file_arch = "../results/Data_out/occs_archetypes_dynamic.csv"
shock_name = "employment_effects_95% by 2035_"


# state whether shock is net or relative to reference scenario
shock_type = "relative"#'net'
# depending on 
if shock_type == "net":
    path_data_scenarios = "../results/Data_out/results_sensitivity/annual_change_sens/"
    # shock_name = "employment_effects_95% by 2035_"
    # file_shock = "employment_effects_95% by 2035_2021-2022__3_oct_2022.csv"
    shock_name = '4_empl_95% by 2035_'
    file_shock = "4_empl_95% by 2035_2020-2021.csv"
    file_out_name = "occ_archetypes_thresholds_net_"
    shock_verison = ".csv"
elif shock_type == "relative":
    path_data_scenarios = "../results/Data_out/results/"
    shock_name = "diff_baseline_employment_effects_"
    file_shock = "diff_baseline_employment_effects_2044-2045__3_oct_2022.csv"
    file_out_name = "occ_archetypes_thresholds_relbase_"
    shock_verison = "__3_oct_2022.csv"
else:
    print("error incorrect shock type")


start_year = 2020
mark_year = 2034 # year marking the period split. Exclusive of mark year
# for first period, inclusive for latter period. 
end_year = 2038
print(mark_year, end_year)

# shock variables
# dividing into opex and capex
technologies_names = ["Biomass", "Coal", "Solar", "Natural gas", "Geothermal", "Hydro", \
                      "Nuclear", "Wind", "Battery storage", "Transmission and distribution"]
all_opex = [name + "-opex" for name in technologies_names]
all_capex = [name + "-capex" for name in technologies_names]
all_capex.remove("Nuclear-capex")
all_capex.remove("Transmission and distribution-capex")
all_capex.append('T&D-capex')

# Get O*NET green skills
df_onet = pd.read_csv(path_data_maria + "Green_Occupations.csv")
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

###############
# Importing shocks
###############
# Import an initial shock file to help define colums
df_shock = pd.read_csv(path_data_scenarios + file_shock)


###############
# Splitting shocks into positive, negative, opex and capex
###############
# Defining new dataframe to put shocks
df_all_shocks = pd.DataFrame()
df_all_shocks['O*NET-SOC Code'] = df_shock['Unnamed: 0']
df_all_shocks['OCC_TITLE'] = df_shock['OCC_TITLE']
df_all_shocks["TOT_EMP"] = df_shock["TOT_EMP"]
df_all_shocks["A_MEAN"] = df_shock["A_MEAN"]
list_soc_occ = list(df_all_shocks['O*NET-SOC Code'])


# Get an estimate of the classification system of shocks
print("Number of occupations at 6 digit", len(df_all_shocks['O*NET-SOC Code']\
    [~df_all_shocks['O*NET-SOC Code'].str.endswith("0")]))
print("Number of occuations at 5 or digits",len(df_all_shocks['O*NET-SOC Code']\
    [df_all_shocks['O*NET-SOC Code'].str.endswith("0")]))
print("Number of occupations at 3 digit",len(df_all_shocks['O*NET-SOC Code']\
    [df_all_shocks['O*NET-SOC Code'].str.endswith("000")]))

# set columns to fill in shocks
# total shocks
df_all_shocks['total_shock'] = 0
df_all_shocks['total_shock_positive'] = 0
df_all_shocks['total_shock_negative'] = 0
# shocks in Opex and Capex
df_all_shocks['total_shock_opex'] =0
df_all_shocks['total_shock_opex_positive'] =0
df_all_shocks['total_shock_opex_negative'] =0
df_all_shocks['total_shock_capex'] =0
df_all_shocks['total_shock_capex_positive'] =0
df_all_shocks['total_shock_capex_negative'] =0

for i in range(start_year, end_year):
    # read df with shock for each year
    df_temp = pd.read_csv(path_data_scenarios + shock_name + \
                          str(i) + "-" + str(i+1) + shock_verison)
    df_temp["shock_opex"] = df_temp[all_opex].sum(axis=1)
    df_temp["shock_capex"] = df_temp[all_capex].sum(axis=1)
    
    df_all_shocks['total_shock'] +=  df_temp ["emp_tot"]
    df_all_shocks['total_shock_opex'] +=  df_temp ["shock_opex"]
    df_all_shocks['total_shock_capex'] +=  df_temp ["shock_capex"]
    # doesn't matter, no overlap it seems
    df_all_shocks['total_shock_positive'] += np.maximum( df_temp ["emp_tot"], 0)
    df_all_shocks['total_shock_negative'] += np.minimum( df_temp ["emp_tot"], 0)
    df_all_shocks['total_shock_opex_positive'] += np.maximum( df_temp ["shock_opex"], 0)
    df_all_shocks['total_shock_opex_negative'] += np.minimum( df_temp ["shock_opex"], 0)
    df_all_shocks['total_shock_capex_positive'] += np.maximum( df_temp ["shock_capex"], 0)
    df_all_shocks['total_shock_capex_negative'] += np.minimum( df_temp ["shock_capex"], 0)
     
    
###############
# Splitting shocks into before and after mark
###############
df_shock_dynamics = pd.DataFrame()
df_shock_dynamics['O*NET-SOC Code'] = df_shock['Unnamed: 0']
df_shock_dynamics['OCC_TITLE'] = df_shock['OCC_TITLE']
df_shock_dynamics["TOT_EMP"] = df_shock["TOT_EMP"]

for i in range(start_year, end_year):
    # read df with shock for each year
    df_temp = pd.read_csv(path_data_scenarios +shock_name + \
                          str(i) + "-" + str(i+1) + shock_verison)

    df_shock_dynamics['shock ' + str(i)] =  df_temp["emp_tot"]
    df_shock_dynamics['shock_opex ' + str(i)] =  df_temp[all_opex].sum(axis=1)
    df_shock_dynamics['shock_capex ' + str(i)] =  df_temp[all_capex].sum(axis=1)
    
# splitting between capex and opex
df_shock_dynamics_capex = pd.DataFrame(columns=list(df_shock['Unnamed: 0']), \
    index=[i for i in range(start_year, end_year)])
df_shock_dynamics_opex = pd.DataFrame(columns=list(df_shock['Unnamed: 0']), \
    index=[i for i in range(start_year, end_year)])

for i in range(start_year, end_year):
    df_temp = pd.read_csv(path_data_scenarios +shock_name + \
                          str(i) + "-" + str(i+1) + shock_verison)
    df_shock_dynamics_capex.loc[i] =  np.array(df_temp[all_capex].sum(axis=1))
    df_shock_dynamics_opex.loc[i] =  np.array(df_temp[all_opex].sum(axis=1))    
    
    
# get shocks before and after for both opex and capex (sum separeteyly)
df_all_shocks["shock_after_" + str(mark_year)] = np.array(df_shock_dynamics_opex\
    .loc[mark_year:].sum(axis=0) + \
    df_shock_dynamics_capex.loc[mark_year:].sum(axis=0))
df_all_shocks["shock_before_" + str(mark_year)] = np.array(df_shock_dynamics_opex\
    .loc[:mark_year].sum(axis=0) + \
    df_shock_dynamics_capex.loc[:mark_year].sum(axis=0))

# Get average positive shock before and after
av_af_pos = df_all_shocks["shock_after_" + str(mark_year)].loc[\
    df_all_shocks["shock_after_" + str(mark_year)] > 0].mean()
av_af_neg = df_all_shocks["shock_after_" + str(mark_year)].loc[\
    df_all_shocks["shock_after_" + str(mark_year)] < 0].mean()
av_bf_neg = df_all_shocks["shock_before_" + str(mark_year)].loc[\
    df_all_shocks["shock_before_" + str(mark_year)] < 0].mean()
av_bf_pos = df_all_shocks["shock_before_" + str(mark_year)].loc[\
    df_all_shocks["shock_before_" + str(mark_year)] > 0].mean()


df_arch_joris = pd.read_csv(file_arch)

df_arch_joris['temp'] = np.where(df_arch_joris['temp'] > 0.01, \
    df_arch_joris['temp'], 0)
df_arch_joris['perm'] = np.where(df_arch_joris['perm'] > 0.01, \
    df_arch_joris['perm'], 0)
df_arch_joris['lost'] = np.where(df_arch_joris['lost'] > 0.01, \
    df_arch_joris['lost'], 0)


df_arch_joris['temp_binary'] = df_arch_joris.apply(lambda x: 1 if (x['temp'] ==\
    max(x['temp'], x['lost'], x['perm']) and x['temp'] > 0.01) else 0, axis=1)
df_arch_joris['lost_binary'] = df_arch_joris.apply(lambda x: 1 if (x['lost'] == \
    max(x['temp'], x['lost'], x['perm']) and x['lost'] > 0.01) else 0, axis=1)
df_arch_joris['perm_binary'] = df_arch_joris.apply(lambda x: 1 if (x['perm'] == \
    max(x['temp'], x['lost'], x['perm']) and  x['perm'] > 0.01) else 0, axis=1)

print("dynamic archaetypes" )
for t in ['temp_binary', 'lost_binary', 'perm_binary']:
    print("occ in ",t, df_arch_joris[t].sum())

def add_columns_archetypes_employment(df, p=0.01):
    '''Classifies occupations into block splitting by before/after shock
    being larger than p*employment with corresponding sign
    '''
    df['Phase_out_p' + str(p)] = np.where((df["shock_before_"+ str(mark_year)] \
        < -p*df["TOT_EMP"]) & (df["shock_after_" + str(mark_year)] \
        < -p*df["TOT_EMP"]), 1, 0)
    df['Permanent_boost_p' + str(p)] = np.where((df["shock_before_"+ str(mark_year)] \
        > p*df["TOT_EMP"]) & (df["shock_after_" + str(mark_year)] \
        > p*df["TOT_EMP"]), 1, 0)
    df['Temporary_boost_p' + str(p)] = np.where((df["shock_before_"+ str(mark_year)] \
        > p*df["TOT_EMP"]) & (df["shock_after_" + str(mark_year)] \
        < -p*df["TOT_EMP"]), 1, 0)
    df['Late_boost_p' + str(p)] = np.where((df["shock_before_"+ str(mark_year)] \
        < -p*df["TOT_EMP"]) & (df["shock_after_" + str(mark_year)] \
        > p*df["TOT_EMP"]), 1, 0)
   
   
def add_columns_archetypes_employment_sqrt(df, r=0.01):
      df['Phase_out_r' + str(r)] = np.where((df["shock_before_"+ str(mark_year)] \
        < 0) & (df["shock_after_" + str(mark_year)] < 0) \
        & ((df["shock_before_"+ str(mark_year)]/df["TOT_EMP"]) ** 2 \
        + (df["shock_after_" + str(mark_year)]/df["TOT_EMP"]) ** 2 \
        > r**2) , 1, 0)  
      df['Permanent_boost_r' + str(r)] = np.where((df["shock_before_"+ str(mark_year)] \
        > 0) & (df["shock_after_" + str(mark_year)] > 0) \
        & ((df["shock_before_"+ str(mark_year)]/df["TOT_EMP"]) ** 2 \
        + (df["shock_after_" + str(mark_year)]/df["TOT_EMP"]) ** 2 \
        > r**2) , 1, 0) 
      df['Temporary_boost_r' + str(r)] = np.where((df["shock_before_"+ str(mark_year)] \
        > 0) & (df["shock_after_" + str(mark_year)] < 0) \
        & ((df["shock_before_"+ str(mark_year)]/df["TOT_EMP"]) ** 2 \
        + (df["shock_after_" + str(mark_year)]/df["TOT_EMP"]) ** 2 \
        > r**2) , 1, 0) 
      df['Late_boost_r' + str(r)] = np.where((df["shock_before_"+ str(mark_year)] \
        < 0) & (df["shock_after_" + str(mark_year)] > 0) \
        & ((df["shock_before_"+ str(mark_year)]/df["TOT_EMP"]) ** 2 \
        + (df["shock_after_" + str(mark_year)]/df["TOT_EMP"]) ** 2 \
        > r**2) , 1, 0) 
   
add_columns_archetypes_employment(df_all_shocks, p = 0.01)  
add_columns_archetypes_employment(df_all_shocks, p = 0.001)
add_columns_archetypes_employment(df_all_shocks, p = 0.005)
add_columns_archetypes_employment_sqrt(df_all_shocks, r = 0.01)
add_columns_archetypes_employment_sqrt(df_all_shocks, r = 0.005)
add_columns_archetypes_employment_sqrt(df_all_shocks, r = 0.001)

archetypes = ['Phase_out_r0.01', 'Permanent_boost_r0.01',\
    'Temporary_boost_r0.01', 'Late_boost_r0.01']

for t in archetypes:
    print("occ in ",t, df_all_shocks[t].sum())

df_all_shocks.to_csv(path_data_maria + file_out_name\
    + str(mark_year) + "_" + str(end_year) + ".csv", \
    index=False)

print('exporting archetypes done')

### Testing results with OPEX and CAPEX differently
df_all_shocks[["shock_before_"+ str(mark_year), "shock_after_" + str(mark_year)]]\
    .loc[df_all_shocks['OCC_TITLE'] == 'Helpers--Extraction Workers']

list_shocks_str_bf  =['shock ' + str(i) for i in range(start_year, mark_year)]
list_shocks_str_af  =['shock ' + str(i) for i in range(mark_year, end_year)]

df_subset = df_shock_dynamics[['OCC_TITLE', 'O*NET-SOC Code'] + list_shocks_str_bf]
df_subset_af = df_shock_dynamics[['OCC_TITLE', 'O*NET-SOC Code'] + list_shocks_str_af]

df_subset.loc[df_subset['OCC_TITLE'] == 'Helpers--Extraction Workers'][list_shocks_str_bf]

df_subset_af.loc[df_subset_af['OCC_TITLE'] == 'Helpers--Extraction Workers'][list_shocks_str_af]


ext_workers_soc = '47-5080'

df_shock_dynamics_opex[ext_workers_soc]
df_shock_dynamics_capex[ext_workers_soc]

archaetypes = ['Phase_out', 'Permanent_boost', 'Temporary_boost', 'Late_boost']

for p in [0.001, 0.005]:
    print("percentage of employment p = ", p )
    for t in archaetypes:
        print("occ in ",t, df_all_shocks[t + "_p" + str(p)].sum())
        
for r in [ 0.01, 0.005,0.001]:
    print("percentage of employment r = ", r )
    for t in archaetypes:
        print("occ in ",t, df_all_shocks[t + "_r" + str(r)].sum())
        
#### Comparing presicion and recall

for r in [ 0.01, 0.005, 0.001]:
    print("sum radious of employment r = ", r )
    print("Phase out precision ", precision_score(df_all_shocks['Phase_out_r' + str(r)], \
        df_arch_joris['lost_binary']), " recall ", \
            recall_score(df_all_shocks['Phase_out_r' + str(r)], \
            df_arch_joris['lost_binary'])    )
    print("Perm boost precision ", precision_score(df_all_shocks['Permanent_boost_r' + str(r)], \
        df_arch_joris['perm_binary']), " recall ", \
            recall_score(df_all_shocks['Permanent_boost_r' + str(r)], \
            df_arch_joris['perm_binary']))
    print("Temp boost precision ", precision_score(df_all_shocks['Temporary_boost_r' + str(r)], \
        df_arch_joris['temp_binary']), " recall ", \
            recall_score(df_all_shocks['Temporary_boost_r' + str(r)], \
            df_arch_joris['temp_binary']))
    for t in archaetypes:
        print("occ in ",t, df_all_shocks[t + "_r" + str(r)].sum())


for p in [ 0.01, 0.005, 0.001]:
    print("percentage of employment p = ", p  )
    print("Phase out precision ", precision_score(df_all_shocks['Phase_out_p' + str(p)], \
        df_arch_joris['lost_binary']), " recall ", \
            recall_score(df_all_shocks['Phase_out_r' + str(r)], \
            df_arch_joris['lost_binary'])    )
    print("Perm boost precision ", precision_score(df_all_shocks['Permanent_boost_p' + str(p)], \
        df_arch_joris['perm_binary']), " recall ", \
            recall_score(df_all_shocks['Permanent_boost_p' + str(p)], \
            df_arch_joris['perm_binary']))
    print("Temp boost precision ", precision_score(df_all_shocks['Temporary_boost_p' + str(p)], \
        df_arch_joris['temp_binary']), " recall ", \
            recall_score(df_all_shocks['Temporary_boost_p' + str(p)], \
            df_arch_joris['temp_binary']))
    for t in archaetypes:
        print("occ in ",t, df_all_shocks[t + "_p" + str(p)].sum())


# dividing occupations into archaetypes
# defining by positive and negative




# restriction to shock being above 1% employment
# both parts
p = 0.01

df_all_shocks['Phase_out'] = np.where((df_all_shocks["shock_before_"\
    + str(mark_year)] < -p*df_all_shocks['TOT_EMP']) & \
        (df_all_shocks["shock_after_" + str(mark_year)] < -p*df_all_shocks['TOT_EMP']), 1, 0)
df_all_shocks['Permanent_boost'] = np.where((df_all_shocks["shock_before_"\
    + str(mark_year)] > p*df_all_shocks['TOT_EMP']) & \
        (df_all_shocks["shock_after_" + str(mark_year)] > p*df_all_shocks['TOT_EMP']), 1, 0)
df_all_shocks['Temporary_boost'] = np.where((df_all_shocks["shock_before_"\
    + str(mark_year)] > p*df_all_shocks['TOT_EMP']) & \
        (df_all_shocks["shock_after_" + str(mark_year)] < -p*df_all_shocks['TOT_EMP']), 1, 0)
df_all_shocks['Late_boost'] = np.where((df_all_shocks["shock_before_"\
    + str(mark_year)] < -p*df_all_shocks['TOT_EMP']) & \
        (df_all_shocks["shock_after_" + str(mark_year)] > p*df_all_shocks['TOT_EMP']), 1, 0)


print("percentage of employment p = ", p )
for t in archaetypes:
    print("occ in ",t, df_all_shocks[t].sum())


    
# restriction to shock being above 1% employment
# both parts
p = 0.001

df_all_shocks['Phase_out'] = np.where((df_all_shocks["shock_before_"\
    + str(mark_year)] < -p*df_all_shocks['TOT_EMP']) & \
        (df_all_shocks["shock_after_" + str(mark_year)] < -p*df_all_shocks['TOT_EMP']), 1, 0)
df_all_shocks['Permanent_boost'] = np.where((df_all_shocks["shock_before_"\
    + str(mark_year)] > p*df_all_shocks['TOT_EMP']) & \
        (df_all_shocks["shock_after_" + str(mark_year)] > p*df_all_shocks['TOT_EMP']), 1, 0)
df_all_shocks['Temporary_boost'] = np.where((df_all_shocks["shock_before_"\
    + str(mark_year)] > p*df_all_shocks['TOT_EMP']) & \
        (df_all_shocks["shock_after_" + str(mark_year)] < -p*df_all_shocks['TOT_EMP']), 1, 0)
df_all_shocks['Late_boost'] = np.where((df_all_shocks["shock_before_"\
    + str(mark_year)] < -p*df_all_shocks['TOT_EMP']) & \
        (df_all_shocks["shock_after_" + str(mark_year)] > p*df_all_shocks['TOT_EMP']), 1, 0)


print("percentage of employment p = ", p )
for t in archaetypes:
    print("occ in ",t, df_all_shocks[t].sum())
# restriction to shock being above mean before/after pos/neg respectively

df_all_shocks['Phase_out'] = np.where((df_all_shocks["shock_before_"\
    + str(mark_year)] < av_bf_neg) & \
        (df_all_shocks["shock_after_" + str(mark_year)] < av_af_neg), 1, 0)
df_all_shocks['Permanent_boost'] = np.where((df_all_shocks["shock_before_"\
    + str(mark_year)] > av_bf_pos) & \
        (df_all_shocks["shock_after_" + str(mark_year)] > av_af_neg), 1, 0)
df_all_shocks['Temporary_boost'] = np.where((df_all_shocks["shock_before_"\
    + str(mark_year)] > av_bf_pos) & \
        (df_all_shocks["shock_after_" + str(mark_year)] < av_af_neg), 1, 0)
df_all_shocks['Late_boost'] = np.where((df_all_shocks["shock_before_"\
    + str(mark_year)] < av_bf_neg) & \
        (df_all_shocks["shock_after_" + str(mark_year)] > av_af_pos), 1, 0)

print("above mean before/after pos/neg")
for t in archaetypes:
    print("occ in ",t, df_all_shocks[t].sum())

########
# Compare with previous 
########

# restriction to shock being above 1% employment
# both parts
p = 0.001

df_all_shocks['Phase_out'] = np.where((df_all_shocks["shock_before_"\
    + str(mark_year)] < -p*df_all_shocks['TOT_EMP']) & \
        (df_all_shocks["shock_after_" + str(mark_year)] < -p*df_all_shocks['TOT_EMP']), 1, 0)
df_all_shocks['Permanent_boost'] = np.where((df_all_shocks["shock_before_"\
    + str(mark_year)] > p*df_all_shocks['TOT_EMP']) & \
        (df_all_shocks["shock_after_" + str(mark_year)] > p*df_all_shocks['TOT_EMP']), 1, 0)
df_all_shocks['Temporary_boost'] = np.where((df_all_shocks["shock_before_"\
    + str(mark_year)] > p*df_all_shocks['TOT_EMP']) & \
        (df_all_shocks["shock_after_" + str(mark_year)] < -p*df_all_shocks['TOT_EMP']), 1, 0)
df_all_shocks['Late_boost'] = np.where((df_all_shocks["shock_before_"\
    + str(mark_year)] < -p*df_all_shocks['TOT_EMP']) & \
        (df_all_shocks["shock_after_" + str(mark_year)] > p*df_all_shocks['TOT_EMP']), 1, 0)


print("percentage of employment p = ", p )
for t in archaetypes:
    print("occ in ",t, df_all_shocks[t].sum())
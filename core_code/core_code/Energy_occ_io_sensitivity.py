import pandas as pd
import numpy as np
global MODELRUN_ID #define here so that it can be global and set below if empty

def run_model(MODELRUN_ID):

    # load modelrun settings    
    import settings
    run_settings = settings.run_settings(MODELRUN_ID)
        
    #save the config settings for this specific model run in the data_out folder
    run_settings.config_file.to_csv(run_settings.FL_data_out + str(MODELRUN_ID) + "_config.csv", index=False)

    import help_sens # includes sensitivity analysis functions to update cost vectors, and T&D calculation
    
    # accomodate the sensitivity analysis for OPEX I/O matrix
    if run_settings.SENS_ANALYSIS_OPEX_COST_VECTORS == True:
        #get an array of sub-runs for sensitivity analysis involving sub_runs
        subrun_array = np.random.randint(low=0, high=1000, size = run_settings.SENS_ANALYSIS_SUBRUNS)
    else:
    # if a seed exists use it (enables a replication of a model run that uses a random number generator)
    # put whatever seed value is used into an array so it can be iterated in loop below
        if run_settings.USE_SEED == True:
            subrun_array = np.full(1, run_settings.SEED)
        else: #otherwise just do one model run with No Seed
            subrun_array = np.full(1, 0)

    for i, sub_run_num in enumerate(subrun_array):
        # set the SEED to the sub_run_num in settings as it is used by IO_scramble in IO_2018_include_opex_IO
        # (so the scrambled run can be recreated if required)
        run_settings.SEED = sub_run_num 
        import IO_2018_include_opex_IO # computes new IO table with opex cost vector sensitivity changes
        IO_assert_error = IO_2018_include_opex_IO.run(MODELRUN_ID, run_settings)
        while IO_assert_error is not None:
            # seed creates values that cause error in merging opex with IO table (e.g. 578). 
            # We generate a new seed and store that
            sub_run_num = np.random.randint(low=0, high=1000)
            subrun_array[i] = sub_run_num
            run_settings.SEED = sub_run_num
            IO_assert_error = IO_2018_include_opex_IO.run(MODELRUN_ID, run_settings)
    
        import TransDist_opexcapex # for sensitivity analysis on T&D calc
        TransDist_opexcapex.run(MODELRUN_ID, run_settings)

        # run the model
        import run_model as run_model
        run_model.run(MODELRUN_ID, run_settings)

    if run_settings.SENS_ANALYSIS_OPEX_COST_VECTORS == True:
        #save the sub_runs seed numbers in the data_out folder
        df_subrun_array = pd.DataFrame(subrun_array)
        df_subrun_array.to_csv(run_settings.FL_data_out + str(MODELRUN_ID) + "_opex_subruns.csv", index=False)

    # Only create figures if the setting is true (and only for the last model sub_run)
    if run_settings.CREATE_FIGURES == True:
        import create_figures_aggregate_jobs
        create_figures_aggregate_jobs.create_figs(MODELRUN_ID, run_settings)

#Load the modelrun file 
## (I have seperated the modelrun_id from the rest of the config file so that I can put it in the .gitignore)
mrid_file = pd.read_csv("modelrun_id.csv")

#  get last model_run_id from modelrun file
# iterate it by one and save the file again with the updaed number
mrid = mrid_file.loc[0,'modelrun_id'].astype(int)
mrid = mrid + 1 #add one to model_run_id so it uses a unique one for each model run
mrid_file.loc[0,'modelrun_id']=mrid #modify the dataframe so the new modelrun_id can be saved

##save the latest model_run_id in the master_config file so the next one will be iterated
mrid_file.to_csv("modelrun_id.csv", index=False) 

MODELRUN_ID = mrid #set global modelrun_id

    
#Load the config file 
config_file = pd.read_csv("master_config.csv")

# run the model with the new modelrun_id
run_model(MODELRUN_ID) 

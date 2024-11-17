import pandas as pd

class run_settings:


    def __init__(self,MODELRUN_ID):
        #set global variables
        self.name_ending = "__3_oct_2022"

        # set global file locations
        self.FL_data_sens = '../../data/data_sens/'
        self.FL_config =""
        self.FL_data  = '../../data/'
        self.FL_energy = '../../energy_scenario_code/data_out/'
        self.FL_analysis  = "../../analysis/"
        self.FL_preprocessed_data = "../../results/data_out/"
        self.FL_ann_data_out = "../../results/data_out/results_sensitivity/annual_change_sens/" 
        self.FL_data_out = "../../results/data_out/results_sensitivity/" #seperate the annual data from the summary data 
        self.FL_fig_folder = "../../results/figs/"
    

        ##### get global variables from master_config file

        #Load the config file 
        self.config_file = pd.read_csv(self.FL_config + "master_config.csv")

        ####create a new modelrun_id
        # do this only once, not everytime the settings module is imported

        # Sensitivity analysis parameters
        self.SENS_ANALYSIS_COST = self.config_file.loc[0,'SENS_ANALYSIS_COST'].astype(bool) #True
        self.COST_SCENARIO = str(self.config_file.loc[0,'COST_SCENARIO']) # usually 'conservative' but can also be 'moderate', or 'advanced'

        # Sensitibity analysis on Cost vectors
        ## SIGMA determins the spread of the noise
        ## SENS_ANALYSIS_SUBRUNS determines how many sub runs to try

        self.SENS_ANALYSIS_CAPEX_COST_VECTORS = self.config_file.loc[0,'SENS_ANALYSIS_CAPEX_COST_VECTORS'].astype(bool) #True
        self.SENS_ANALYSIS_OPEX_COST_VECTORS = self.config_file.loc[0,'SENS_ANALYSIS_OPEX_COST_VECTORS'].astype(bool) #True
        # you can't run with both true so set OPEX to false
        if self.SENS_ANALYSIS_CAPEX_COST_VECTORS == True:
            self.SENS_ANALYSIS_OPEX_COST_VECTORS = False
                
        self.SIGMA = self.config_file.loc[0,'SIGMA'].astype(float) #0.5
        #number of sub runs for a sensitivity analysis that uses subruns e.g. noise
        self.SENS_ANALYSIS_SUBRUNS= self.config_file.loc[0,'SENS_ANALYSIS_SUBRUNS'].astype(int) #30
        #Determine if a random seed is being used (to replicate a result set with noise) and get that seed
        self.USE_SEED = self.config_file.loc[0,'USE_SEED'].astype(bool) #False
        # if useing the seed then turn remember the seed and off the sensitivity analysis because you would then just get multiple runs using the same seed value
        if self.USE_SEED==True:
            self.SEED = self.config_file.loc[0,'SEED'].astype(int) 
            self.SENS_ANALYSIS_CAPEX_COST_VECTORS = False
            self.SENS_ANALYSIS_OPEX_COST_VECTORS = False
        else: #else set SEED to 0 so code knows not to use a Seed
            self.SEED = 0

        # smoothing
        self.SENS_ANALYSIS_SMOOTHING = self.config_file.loc[0,'SENS_ANALYSIS_SMOOTHING'].astype(bool) #True
        self.YEARS_SMOOTHING_DEFAULT = 5
        self.YEARS_SMOOTHING_SENS = self.config_file.loc[0,'YEARS_SMOOTHING_SENS'].astype(int) # default is 5 but can vary between 1 (no smoothing), 3, 5, and 7?

        # IO year
        self.SENS_ANALYSIS_IO_YEAR = self.config_file.loc[0,'SENS_ANALYSIS_IO_YEAR'].astype(bool) #True
        self.IO_YEAR_DEFAULT = 2018
        self.IO_YEAR_SENS = self.config_file.loc[0,'IO_YEAR_SENS'].astype(int) # default 2018, can vary 2015-2019 incl

        # T&D
        self.SENS_ANALYSIS_TANDD = self.config_file.loc[0,'SENS_ANALYSIS_TANDD'].astype(bool) #True
        # per mile line costs: set default and get sens value
        # usually set to 1433 but vary between 932 - 3624: discounted values following footnote 18 of https://www.nrel.gov/docs/fy17osti/67240.pdf
        self.LINE_COST_DEFAULT = 1433
        self.LINE_COST_SENS = self.config_file.loc[0,'LINE_COST_SENS'].astype(int) #1433

        self.SENS_ANALYSIS_OCC  = self.config_file.loc[0,'SENS_ANALYSIS_OCC'].astype(bool) #True

        #set defaults and get sens values
        self.THRICE_LINE_COST_FACTOR_DEFAULT = 1.37
        self.THRICE_LINE_COST_FACTOR_SENS  = self.config_file.loc[0,'THRICE_LINE_COST_FACTOR_SENS'].astype(float) #1.37
        self.THRICE_LINE_OPEX_FACTOR_DEFAULT  = 1.37
        self.THRICE_LINE_OPEX_FACTOR_SENS = self.config_file.loc[0,'THRICE_LINE_OPEX_FACTOR_SENS'].astype(float) #1.37

        # import export sensitivity
        self.SENS_IMPEXP = self.config_file.loc[0, 'SENS_ANALYSIS_IMPORTEXPORT'].astype(bool) #False
        self.IMPEXP_DEFAULT = 'fixed'
        self.IMPEXP_SENS_OPTION = self.config_file.loc[0, 'IMPORTEXPORT_SENS_OPTION'] # default is 'fixed' but can also be 'decreasing_imports', 'increasing_exports', or 'both', overwrites IMPEXP_DEFAULT

        # T&D construction breakdown sensitivity
        self.SENS_TANDD_CONSTR_BREAKDOWN = self.config_file.loc[0,'SENS_TANDD_CONSTRUCTION'].astype(bool)

        # zero or double solar and wind sensitivity
        self.SENS_ZERO_DOUBLE_SOLAR_WIND = self.config_file.loc[0,'SENS_ZERO_DOUBLE_SOLAR_WIND'].astype(bool)
        self.SOLAR_WIND_SENS_OPTION = self.config_file.loc[0,'SOLAR_WIND_SENS_OPTION']
        
        # Occupation sensitivity
        self.OCC_DEFAULT_OPTION = 'standard'
        self.OCC_SENS_OPTION = str(self.config_file.loc[0,'OCC_SENS_OPTION']) # default is 'standard' but can also be 'low' or 'high' , overwrites OCC_DEFAULT_OPTION

        self.CREATE_FIGURES = self.config_file.loc[0,'CREATE_FIGURES'].astype(bool) #False

        #Set some global variables that are based on config variables
        self.OCC_DEFAULT_PREFIX = 'occ_breakdown/' + self.OCC_DEFAULT_OPTION + '_'
        self.OCC_SENS_PREFIX = 'occ_breakdown/' + self.OCC_SENS_OPTION + '_'
        
        #keep track of the number of sub_runs numbers - iterated in run_model.py
        self.sub_run_id = 0
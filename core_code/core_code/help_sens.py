from audioop import error
import numpy as np
import pandas as pd

def scramble_IO(A, settings):
    '''Multiply all items of matrix A with (1+x) with x~N(0, sigma)
    input: A is a pd.DataFrame 
    input: sigma is a scalar, default sigma at top of this file
    output: B is a pd.Dataframe, with multiplicative noise added, but column sums preserved
    '''


    #set the seed value so we can replicate the results if we want to
    if settings.SEED != None and settings.SEED != 0:
        np.random.seed(settings.SEED)

    # noise is random normal, clipped to be non-zero positive
    noise = np.clip(1 + np.random.normal(0, settings.SIGMA, size=A.shape), 0.001, None)

    B = A.multiply(noise)

    # column sum should be preserved
    B_sum = B.sum(axis=0)
    A_sum = A.sum(axis=0)

    B = B.divide(B_sum, axis=1).multiply(A_sum, axis=1)
    B.fillna(0.0, inplace=True) # if rowsums are zero, div/zero will give n/a: fill with 0s

    return B


def TandD_sens(settings):
    ''' CHECK GLOBAL VARIABLES IN SETTINGS.PY
    
    assign different values to T&D calculation for sensitivity analysis
    output: line_cost = cost per MW-mi of transmission cable (default 1433 USD(2018))
    output: thrice_line_cost_factor = capex cost factor of upgrading lines with 3x capacity (default 1.37)
    output: thrice_line_opex_factor = opex cost factor of upgrading lines with 3x capacity (default 1.37)
    '''
    
    if settings.SENS_ANALYSIS_TANDD == True:
        line_cost,thrice_line_cost_factor, thrice_line_opex_factor = (settings.LINE_COST_SENS,
                                                                  settings.THRICE_LINE_COST_FACTOR_SENS, 
                                                                  settings.THRICE_LINE_OPEX_FACTOR_SENS)
    else:
        line_cost,thrice_line_cost_factor, thrice_line_opex_factor = (settings.LINE_COST_DEFAULT,
                                                                      settings.THRICE_LINE_COST_FACTOR_DEFAULT, 
                                                                      settings.THRICE_LINE_OPEX_FACTOR_DEFAULT)

    return line_cost,thrice_line_cost_factor, thrice_line_opex_factor


def cost_sens(settings):
    ''' CHECK GLOBAL VARIABLE COST_SCENARIO IN SETTINGS.PY
    
    output: file name prefix for right cost trajectory file
    '''
    if settings.COST_SCENARIO == 'advanced':
        file_name_start = 'ADVANCED_'
    elif settings.COST_SCENARIO == 'conservative':
        file_name_start = 'CONSERVATIVE_'
    elif settings.COST_SCENARIO == 'pro-ff':
        file_name_start = 'pro-ff_'
    elif settings.COST_SCENARIO == 'pro-re':
        file_name_start = 'pro-re_'
    elif settings.COST_SCENARIO == 'moderate': # moderate is default
        file_name_start = ''

    return file_name_start


def occ_sens(ind_occ, prse, settings):
    ''' adds or subtracts one std dev from industry-occupation matrix
    
    input: ind_occ = pd dataframe of indusry-occupation employment by BLS
    input: prse = pd dataframe of prse of industry-occupation employment by BLS
    output: ind_occ_sens = pd dataframe to be used for sensitivity analysis
    '''        
    
    if settings.OCC_SENS_OPTION == 'standard':
        ind_occ_sens = ind_occ
    elif settings.OCC_SENS_OPTION == 'high':
        ind_occ_sens = ind_occ.multiply(1 + (prse / 100))
    elif settings.OCC_SENS_OPTION == 'low':
        ind_occ_sens = np.clip(ind_occ.multiply(1 - (prse / 100)), 0, None)
    else: #, catch the error if it is not one of these values
        raise ValueError('settings.OCC_SENS_OPTION has an unknown value') 
    return ind_occ_sens

def switch_year_reduce(year):
    if year < 2030:
        return 1
    elif year < 2035:
        return 0.9
    elif year < 2040:
        return 0.7
    else: # year >= 2040
        return 0.5

def imp_sens(MODELRUN_ID, year, settings):
    ''' reads in correct file given the year and sensitivity analysis
    year: int, the year of the data
    settings: settings object
    '''
    imp_red = switch_year_reduce(year)

    if (settings.IMPEXP_SENS_OPTION == 'fixed') | (settings.IMPEXP_SENS_OPTION == 'increasing_exports') | (imp_red == 1):

        import_perc = pd.read_excel(settings.FL_data + "Data_out_IO/2018-import_perc_sum.xlsx", index_col=0)
        # A_sum_2018_dom = pd.read_csv(
        #     settings.FL_analysis
        #     + "Energy sector disaggregation/Scripts/output_lit/A_2018_inc_elec_split_lit.csv",
        #     index_col=0)
        
        # Z_2018 = pd.read_csv(
        #     settings.FL_analysis
        #     + "Energy sector disaggregation/Scripts/output_lit/Z_2018_inc_elec_split_lit.csv",
        #     index_col=0)


    elif (settings.IMPEXP_SENS_OPTION == 'decreasing_imports') | (settings.IMPEXP_SENS_OPTION == 'both'):

        import_perc = pd.read_excel(settings.FL_data + "Data_out_IO/imp_reduced/2018-import_perc_sum_" + str(imp_red) + ".xlsx", index_col=0)
        # A_sum_2018_dom = pd.read_csv(
        #     settings.FL_data_sens + str(MODELRUN_ID) + '_'
        #     + "A_2018_inc_elec_split_lit_" + str(imp_red) + ".csv",
        #     index_col=0)
        
        # Z_2018 = pd.read_csv(
        #     settings.FL_data_sens + str(MODELRUN_ID) + '_'
        #     + "Z_2018_inc_elec_split_lit_" + str(imp_red) + ".csv",
        #     index_col=0)

    else:
        raise ValueError('settings.IMPEXP_SENS_OPTION has an unknown value')

    return import_perc#, A_sum_2018_dom, Z_2018,

def switch_year_exports(year):
    if year < 2030:
        return 1
    elif year < 2035:
        return 1.1
    elif year < 2040:
        return 1.3
    else: # year >= 2040
        return 1.5
    
def exp_sens(MODELRUN_ID, year, cap_2030_shock, import_p, inv_weights, global_scen, settings):
    
    if (settings.IMPEXP_SENS_OPTION == 'increasing_exports') | (settings.IMPEXP_SENS_OPTION == 'both'):
        exp_inc = switch_year_exports(year)

        export_factor = pd.read_excel(settings.FL_data + 'export_sectors.xlsx', index_col=1)
        # processing
        export_factor.index = export_factor.index.astype(str)
        export_factor.drop(index='22', inplace=True)
        export_factor.drop(index='nan', inplace=True)
        export_factor = export_factor['export']
        # replace Y in 'export' for exp_inc
        export_factor = export_factor.replace('Y', exp_inc-1)
        export_factor = export_factor.replace('N', 1-1)

        inv_weights.drop("22", inplace=True)

        # Solar
        solar_inv_weights = inv_weights[('Solar', 'capex')].multiply(
                    1.0 - import_p.loc[inv_weights.index], axis=0
        )
        solar_inv_weights = solar_inv_weights.multiply(export_factor, axis=0)
        solar_inv_weights *= cap_2030_shock.loc[global_scen, 'Solar'] * 1000

        # Wind
        wind_inv_weights = inv_weights[('Wind', 'capex')].multiply(
                    1.0 - import_p.loc[inv_weights.index], axis=0
        )
        wind_inv_weights = wind_inv_weights.multiply(export_factor, axis=0)
        wind_inv_weights *= cap_2030_shock.loc[global_scen, 'Wind'] * 1000

        return wind_inv_weights, solar_inv_weights
    else:
        return None, None

# def exp_sens(MODELRUN_ID, year, settings):
#     ''' introduces export values for sensitivity analysis
#     year: int, the year of the data
#     settings: settings object
#     '''
#     if (settings.IMPEXP_SENS_OPTION == 'increasing_exports') | (settings.IMPEXP_SENS_OPTION == 'both'):
#         exp_inc = switch_year_exports(year)

#         export_factor = pd.read_excel(settings.FL_data + 'export_sectors.xlsx', index_col=1)
#         # processing
#         export_factor.index = export_factor.index.astype(str)
#         export_factor.drop(index='22', inplace=True)
#         export_factor.drop(index='nan', inplace=True)
#         export_factor = export_factor['export']
#         # replace Y in 'export' for exp_inc
#         export_factor = export_factor.replace('Y', exp_inc-1)
#         export_factor = export_factor.replace('N', 1-1)

#         exports should be independent of the domestic requirements. So we should not 
#         multiply the exports with the domestic requirements, but do it somehow else

#         return export_factor
#     else:
#         return None

def tandd_constr_sens(construction_spending, emp_per_mill, settings):
    ''' introduces construction spending and employment per million for sensitivity analysis
    construction_spending: float, the construction spending
    emp_per_mill: pd dataframe, the employment per million for construction
    settings: settings object
    '''

    # total construction workers per million output
    construction_workers_per_million = emp_per_mill['23'].sum()

    # breakdown of T&D construction occupational employment per workers
    tandd_constr_breakdown = pd.read_csv(settings.FL_preprocessed_data + 'emp_breakdown_TandD_construction_2018.csv', index_col=0)['0']

    # spending * workers / spending * occuaptional employment / worker = occupational employment
    tandd_construction_emp = construction_spending * construction_workers_per_million * tandd_constr_breakdown

    return tandd_construction_emp
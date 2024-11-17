
from ast import Assert
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.testing import assert_array_almost_equal
import help_sens

def run(MODELRUN_ID, settings, imp_reduce = 1):

    # # This notebook disaggregates the 2018 US IO Z matrix utilities sector
    # A. The 2012 US IO Z matrix with detailed industries is used to disaggregate the utilities sector into an electricity sector (221100), a gas sector (221200), and a water sector (221300). 

    # B. A literature based IO matrix is used to disaggregate the electricity sector (221100) further into electric power and control, electric power distribution, and 8 electricity generation technologies (Biomass, fossil fuel, geotherman, hydro, nuclear, solar, wind, and other). The BEA detailed industry output figures for 2018 were used to get the correct split between them.

    # ## A. 2012 detailed BEA IO table interstion into 2018 summary BEA IO table
    # 1. We first read the BEA IO tables
    # 2. Then we take the utility rows and columns to insert and transform them into 71 industry rows and columns
    # 3. Then we make the imputation by tackling the 4 constraints
    # 4. Finally, we combine everything into a 2018 extended IO table

    ### 1. Read in BEA IO tables

    io_year = settings.IO_YEAR_DEFAULT
    if settings.SENS_ANALYSIS_IO_YEAR:
        io_year = settings.IO_YEAR_SENS
    if imp_reduce == 1:
        A_2018 = pd.read_excel(settings.FL_data + "Data_out_IO/" + str(io_year) + "-A_sum_dom.xlsx", index_col=0)
        Z_2018 = pd.read_excel(settings.FL_data + "Data_out_IO/" + str(io_year) + "-Z_sum_dom.xlsx", index_col=0)
        va_2018 = pd.read_excel(settings.FL_data + "Data_out_IO/" + str(io_year) + "-va_sum_dom.xlsx", index_col=0)
        fd_2018 = pd.read_excel(settings.FL_data + "Data_out_IO/" + str(io_year) + "-f_sum_dom.xlsx", index_col=0)

        A_2012 = pd.read_excel(settings.FL_data + "Data_out_IO389/2012-A_det_dom.xlsx", index_col=0)
        Z_2012 = pd.read_excel(settings.FL_data + "Data_out_IO389/2012-Z_det_dom.xlsx", index_col=0)
        va_2012 = pd.read_excel(settings.FL_data + "Data_out_IO389/2012-va_det_dom.xlsx", index_col=0)
        fd_2012 = pd.read_excel(settings.FL_data + "Data_out_IO389/2012-f_det_dom.xlsx", index_col=0)
    else:
        A_2018 = pd.read_excel(settings.FL_data + "Data_out_IO/imp_reduced/" + str(io_year) + "-A_sum_dom_" + str(imp_reduce) + ".xlsx", index_col=0)
        Z_2018 = pd.read_excel(settings.FL_data + "Data_out_IO/imp_reduced/" + str(io_year) + "-Z_sum_dom_" + str(imp_reduce) + ".xlsx", index_col=0)
        va_2018 = pd.read_excel(settings.FL_data + "Data_out_IO/imp_reduced/" + str(io_year) + "-va_sum_dom_" + str(imp_reduce) + ".xlsx", index_col=0)
        fd_2018 = pd.read_excel(settings.FL_data + "Data_out_IO/imp_reduced/" + str(io_year) + "-f_sum_dom_" + str(imp_reduce) + ".xlsx", index_col=0)

        A_2012 = pd.read_excel(settings.FL_data + "Data_out_IO389/imp_reduced/2012-A_det_dom_" + str(imp_reduce) + ".xlsx", index_col=0)
        Z_2012 = pd.read_excel(settings.FL_data + "Data_out_IO389/imp_reduced/2012-Z_det_dom_" + str(imp_reduce) + ".xlsx", index_col=0)
        va_2012 = pd.read_excel(settings.FL_data + "Data_out_IO389/imp_reduced/2012-va_det_dom_" + str(imp_reduce) + ".xlsx", index_col=0)
        fd_2012 = pd.read_excel(settings.FL_data + "Data_out_IO389/imp_reduced/2012-f_det_dom_" + str(imp_reduce) + ".xlsx", index_col=0)

    BEA_71_400 = pd.read_excel(settings.FL_data + "BEA_71_400_naics.xlsx", usecols=[2, 6], index_col=0)
    BEA_71_400.loc['imports'] = 'imports'
    crosswalk_71_400 = BEA_71_400.copy()
    crosswalk_71_400['val'] = 1
    crosswalk_71_400 = crosswalk_71_400.pivot(columns='Detail').fillna(0).droplevel(0, axis=1)

    ### 2. We take the ultily rows and columns of the 2012 400 industry version and crosswalk them into a 71 industry version 

    Z_va_2012 = pd.concat([Z_2012.iloc[:, 21:24], (va_2012.iloc[:, 21:24])])
    Z_va_2018 = pd.concat([Z_2018.iloc[:, [5]], (va_2018.iloc[:, [5]])])

    # groupby the 2012 detailed IO table by summary (71 industry) classification
    Z_va_2012_summary = Z_va_2012.merge(BEA_71_400, left_index=True, right_on='Detail')
    Z_va_2012_summary = Z_va_2012_summary.drop('Detail', axis=1).groupby('Summary').sum()
    Z_va_2012_summary.index = Z_va_2012_summary.index.astype(str)
    Z_va_2012_summary = Z_va_2012_summary.loc[Z_va_2018.index]

    #A_va_2012_summary = Z_va_2012_summary.div(Z_va_2012_summary.sum(axis=0))
    #tot_out_2012 = Z_va_2012.sum()
    #w_2012 = tot_out_2012.div(tot_out_2012.sum())

    #A_va_2018 = Z_va_2018.div(Z_va_2018.sum(axis=0))

    Z_fd_2012 = pd.concat([Z_2012.iloc[21:24, :], fd_2012.iloc[21:24, :]], axis=1)
    Z_fd_2018 = pd.concat([Z_2018.iloc[[5], :], fd_2018.iloc[[5], :]], axis=1)

    Z_fd_2012_summary = Z_fd_2012.T.merge(BEA_71_400, left_index=True, right_on='Detail')
    Z_fd_2012_summary = Z_fd_2012_summary.drop('Detail', axis=1).groupby('Summary').sum().T
    Z_fd_2012_summary.columns = Z_fd_2012_summary.columns.astype(str)
    Z_fd_2012_summary = Z_fd_2012_summary.loc[:, Z_fd_2018.columns]

    Z_core = Z_2012.iloc[21:24, 21:24]

    # ### 3. We will impute the 2012 values into the 2018 IO table. A number of conditions must hold:
    # i. Sum of inputs for disaggregated sectors equals inputs of original

    # ii. Total sum of input from disaggregated sectors equals input from original

    # iii. inputs to and from disaggregated dsectors equals original self-link

    # iV. All items > 0 (including final demand and value added)

    #### i. sum of inputs equals original inputs 

    factor_2018_2012 = Z_va_2018.div(Z_va_2012_summary.sum(axis=1), axis=0)

    # all zeros, so ok!
    zero_list = Z_va_2018[factor_2018_2012.isna()['22']].index
    assert(Z_va_2012_summary.loc[zero_list].sum().sum() ==0)
    assert(Z_va_2018.loc[zero_list].sum().sum() ==0)

    # factor_2018_2012 = A_va_2018.div(A_va_2012_summary.dot(w_2012), axis=0).fillna(0)
    factor_2018_2012 = factor_2018_2012.fillna(0)

    Z_column_insert = Z_va_2012_summary.multiply(factor_2018_2012['22'], axis=0)

    #### ii. sum of outputs = original output

    factor_2018_2012 = Z_fd_2018.div(Z_fd_2012_summary.sum(axis=0), axis=1)

    # all zeros, so no problem
    zero_list = Z_fd_2018.T[factor_2018_2012.isna().loc['22']].T.columns
    assert(Z_fd_2012_summary.loc[:, zero_list].sum().sum() == 0)
    assert(Z_fd_2018.loc[:, zero_list].sum().sum() == 0)

    # factor_2018_2012 = A_va_2018.div(A_va_2012_summary.dot(w_2012), axis=0).fillna(0)
    factor_2018_2012 = factor_2018_2012.fillna(0)

    Z_row_insert = Z_fd_2012_summary.multiply(factor_2018_2012.loc['22'], axis=1)

    # Z_row_insert



    #### iii. Self links consistent


    # ipfp
    def ipfp(b, u, v, X, rowsum=True):
    
        a_n = u / X.dot(b)
        b_n = v / X.T.dot(a_n)
        if rowsum:
            a_n = u / X.dot(b_n)
    
        m_n = np.outer(a_n, b_n)
        m_n = np.multiply(m_n, X)
    
        return a_n, b_n, m_n

    v = Z_column_insert.loc['22'].to_numpy()
    u = Z_row_insert['22'].to_numpy()
    X = Z_core.to_numpy()
    b = np.array([1, 1, 1])

    a_n, b_n, m_n = ipfp(b, u, v, X, rowsum=False)

    np.testing.assert_allclose(m_n.sum(axis=0), v)
    np.testing.assert_allclose(m_n.sum(axis=1), u)

    Z_core_insert = pd.DataFrame(m_n, columns=Z_core.columns, index=Z_core.index)



    #### iv. All positive


    ### 4. Combine everything
    Z_2018_new = pd.concat([
        pd.concat([Z_2018.iloc[:5, :5], Z_column_insert.iloc[:5, :], Z_2018.iloc[:5, 6:]], axis=1), 
        pd.concat([Z_row_insert.iloc[:, :5], Z_core_insert, Z_row_insert.iloc[:, 6:71]], axis=1), 
        pd.concat([Z_2018.iloc[6:, :5], Z_column_insert.iloc[6:72, :], Z_2018.iloc[6:, 6:]], axis=1)],
        axis=0)

    fd_2018_new = pd.concat([fd_2018.iloc[:5, :], Z_row_insert.iloc[:, 71:], fd_2018.iloc[6:, :]], axis=0)

    va_2018_new = pd.concat([va_2018.iloc[:, :5], Z_column_insert.iloc[72:, :], va_2018.iloc[:, 6:]], axis=1)

    Z_2018_new.index = Z_2018_new.index.astype(str)
    va_2018_new.index = va_2018_new.index.astype(str)
    fd_2018_new.index = fd_2018_new.index.astype(str)

    Z_2018_new.columns = Z_2018_new.columns.astype(str)
    va_2018_new.columns = va_2018_new.columns.astype(str)
    fd_2018_new.columns = fd_2018_new.columns.astype(str)

    Z_fd_2018_new = pd.concat([Z_2018_new, fd_2018_new], axis=1)

    Z_va_2018_new = pd.concat([Z_2018_new, va_2018_new])

    Z_va_fd_2018_new = pd.concat([Z_fd_2018_new, va_2018_new])

    #Z_va_fd_2018_new.to_csv('output_lit/Z_2018_incl_utilities_split.csv')

    A_va_2018_new = Z_va_2018_new / Z_va_2018_new.sum()
    #A_va_2018_new.to_csv('output_lit/A_2018_incl_utilities_split.csv')



    # ## B: We next use literature-based estimates to further disaggregate the electricity sector
    # electricity sector (221100)
    # 1. We first read the li IO tables and find the US submatrix
    # 2. We make a crosswalk from EXIOBASE industries to BEA industries
    # 3. We prepare the Z matrix electricity rows and columns from EXIOBASE data using the crosswalk
    # 4. Then we make the imputation by tackling the 4 constraints
    # 5. Finally, we combine everything into a 2018 extended IO table

    ### 1. Read in literature estimates tables and find US submatrix

    A_lit= pd.read_excel(settings.FL_data + 'Synthetic IO sectors energy technologies alt.xlsx', \
                    sheet_name='opex_for_A_inclusion', index_col=1)
    A_lit.drop('Industries', axis=1, inplace=True)
    # Monetary IO 2011
    #A = pd.read_csv('../Data/Exiobase/IOT_2011_ixi/IOT_2011_ixi/A.txt', delimiter = "\t", index_col=[0,1], header=[0,1])
    #fc = pd.read_csv('../Data/Exiobase/IOT_2011_ixi/IOT_2011_ixi/Y.txt', delimiter = "\t", index_col=[0,1], header=[0,1])

    ind_names = {'Wind': 'Wind electric power generation', 'Solar': 'Solar electric power generation', \
                 'Fossil Fuel': 'Fossil fuel electric power generation', 'Biomass': 'Biomass electric power generation', \
                 'Geothermal': 'Geothermal electric power generation', 'Hydro': 'Hydroelectric power generation', \
                 'Coal': 'Coal electric power generation', 'Gas': 'Gas electric power generation'}

    A_lit.rename(ind_names, axis=1, inplace=True)
    A_lit.index = A_lit.index.astype(str)
    A_lit.drop('Fossil fuel electric power generation', axis=1, inplace=True)

    # odd typo and cannot get rid off:
    A_lit.rename({'521CI ': '521CI'}, axis=0, inplace=True)

    # manual edits (expert opinion):
    # Add 2% of opex spending in coal/gas to utilities (currently 0):
    A_lit.loc[:, 'Coal electric power generation'] *= 0.98
    A_lit.loc[:, 'Gas electric power generation'] *= 0.98
    A_lit.loc['22', 'Coal electric power generation'] = 0.02
    A_lit.loc['22', 'Gas electric power generation'] = 0.02
    A_lit_s = A_lit.copy()


    # for those we have no litearture estimates, we impute the values of the electricity sector (221100)
    A_lit_others = pd.DataFrame()
    not_lit = ['Nuclear electric power generation', 'Other electric power generation', \
               'Electric bulk power transmission and control', 'Electric power distribution']

    for sect in not_lit:
        A_lit_others[sect] = A_va_2018_new['221100']

    # except for fuel costs, which we assume to be 1/5th for nuclear (mining), and 0 for the others:
    # nuclear also has fuel costs, but these are not part of the I-O table for two reasons
    # 1: The 2012 I-O table shows no relation between elec geneartion and sector 2122A0 (gold and other 
    # metal mining, incl uranium). 2. The EIA suggests most uranium is imported and about 1/6th is 
    # sourced locall https://www.eia.gov/totalenergy/data/monthly/pdf/sec8_5.pdf
    # Total costs of nuclear fuel for 2018 seems to be 42.98 * 11.1 = 477 million USD (about 1.3% for nuclear A matrix)
    # Nuclear fuel is sector 325189 (other inorganic chemicals)
    nuclear_fuel = 0.013
    for sect in not_lit:
        for fuel_sect in ['211', '212', '324', '486']:
            A_lit_others.loc[fuel_sect, sect] = 0
    A_lit_others.loc['325', 'Nuclear electric power generation'] = nuclear_fuel

    for sect in not_lit:
        A_lit_others.iloc[:74, :].loc[:, sect] /= (A_lit_others.iloc[:74, :].loc[:, sect].sum() / \
                                                (1-A_lit_others.loc[['V001', 'V002', 'V003'], sect].sum()))

    # we lack value added fractions of total spending for the new sectors
    # we impute values of the electricity sector (221100) for v002 and v003
    va_frac = ((va_2018_new.sum() / (va_2018_new.sum() + Z_2018_new.sum())).loc['221100'])

    #A_va_2018_new.loc['212', '221100']

    A_lit = A_lit_others.join(A_lit * (1-va_frac))    



    # similarlty, where we do have data, we also do not have the value added split
    # We impute that with data from the 'other electric power generation'
    ref = 'Other electric power generation'
    for col in A_lit.columns:
        for v in ['V001', 'V002', 'V003']:
            A_lit.loc[v, col] = A_lit.loc[v, ref]

    # The original inputs from the 'utility sector' (22) are replaced with that of the disaggregated sectors)
    elec_sects = ['221100', '221200', '221300']
    elec_sum = A_lit.loc[elec_sects, ref].sum()
    for col in A_lit_s.columns:
        for s in elec_sects:
            A_lit.loc[s, col] = A_lit_s.loc['22', col] * (A_lit.loc[s, ref] / elec_sum) * (1-va_frac)
    

    A_lit = A_lit.fillna(0)

    #set(A_lit_s.index) - set(A_lit.index)



    # BEA total output in 2018 of detailed electricty sectors, in millions of dollar
    bea_electricity_output = pd.read_excel("../../data/BEA/GrossOutput.xlsx", sheet_name='UGO305-A', usecols=[1, 21,22,23,24,25], skiprows=list(range(7)) + list(range(8, 29)), nrows=10)
    bea_electricity_output.set_index('Industry', inplace=True)
    eia_elec_output = pd.read_csv('../../data/US_generated_electricity.csv', index_col=1)

    # BEA output only has 'fossil fuel', not coal and gas separately. We use
    # EIA data with total GWh output to disaggregate the two, assuming relative monetary output
    # reflects relative GWh output

    coal_frac_fossil = (eia_elec_output.loc[io_year, 'coal']) / \
        (eia_elec_output.loc[io_year, 'oil'] + eia_elec_output.loc[io_year, 'gas'] + eia_elec_output.loc[io_year, 'coal'])

    bea_electricity_output.loc['Coal electric power generation', str(io_year)] = (coal_frac_fossil * \
                    bea_electricity_output.loc['Fossil fuel electric power generation', str(io_year)]).round()
    bea_electricity_output.loc['Gas electric power generation', str(io_year)] = ((1-coal_frac_fossil) * \
                    bea_electricity_output.loc['Fossil fuel electric power generation', str(io_year)]).round()
    bea_electricity_output.drop('Fossil fuel electric power generation', axis=0, inplace=True)

    #coal_frac_fossil

    factor_bea_elec = Z_va_2018_new.sum().loc['221100'] / bea_electricity_output[str(io_year)].sum()
    bea_electricity_output.loc[:, str(io_year)] *= factor_bea_elec
    #print(factor_bea_elec)

    # #### i.b: value added v001 correction based on wages
    # 1. We use total wages to compute the value added fraction
    # 2. We assume that V001 (employee compensation) scales with total wages
    # 3. We readjust intermediate spending to compensate
    # 3. The others V002 (taxes and subsidies) and V003 (gross margin) are assumed to be constant

    # jobs_wages has number of workers and mean wage for electricity sector
    jobs_wages = pd.read_csv(settings.FL_preprocessed_data + 'jobs_per_MWh_through_time.csv')
    jobs_wages = jobs_wages[jobs_wages['year'] == io_year]
    jobs_wages['tot_wage'] = jobs_wages['tot_emp'] * jobs_wages['a_mean'] 

    ind_names_l = {'wind': 'Wind electric power generation', 'solar': 'Solar electric power generation', \
                 'fossil': 'Fossil fuel electric power generation', 'biomass': 'Biomass electric power generation', \
                 'geothermal': 'Geothermal electric power generation', 'hydro': 'Hydroelectric power generation', \
                 'nuclear': 'Nuclear electric power generation', 'other': 'Other electric power generation', \
                 'coal': 'Coal electric power generation', 'gas': 'Gas electric power generation', 'total': 'Total'}

    jobs_wages['name'] = jobs_wages['short name'].map(ind_names_l)
    tot_wage = jobs_wages.set_index('name')['tot_wage']

    # for electric power distribution and transmission we don't have separate estimates
    # We assume everyone working in electricity but not in generation works in transmission and distribution
    tot_wage.loc['Electric power transmission and distribution'] = 2*tot_wage['Total'] - tot_wage.sum()
    tot_wage.drop('Total', inplace=True)

    # We split fossil fuel in coal and gas according to GWh generation split
    tot_wage['Coal electric power generation'] = coal_frac_fossil * tot_wage['Fossil fuel electric power generation']
    tot_wage['Gas electric power generation'] = (1-coal_frac_fossil) * tot_wage['Fossil fuel electric power generation']

    # We further split distribution and transmission by using total USD output data
    frac_dist = bea_electricity_output.loc['Electric power distribution', str(io_year)] / (bea_electricity_output.loc[\
                        'Electric bulk power transmission and control', str(io_year)] + bea_electricity_output.loc[\
                        'Electric power distribution', str(io_year)])

    tot_wage['Electric power distribution'] = frac_dist * tot_wage['Electric power transmission and distribution']
    tot_wage['Electric bulk power transmission and control'] = (1-frac_dist) * tot_wage['Electric power transmission and distribution']

    # we rescale total wage to million of USD, and calculate fraction of total output
    output_vs_wage = pd.concat([bea_electricity_output.loc[:, [str(io_year)]], tot_wage/1000000], axis=1).dropna()
    output_vs_wage['frac_wage_of_output'] = output_vs_wage.tot_wage / output_vs_wage.loc[:,str(io_year)]
    va001_frac_output = (va_2018_new.loc['V001'] / (va_2018_new.sum() + Z_2018_new.sum())).loc['221100']
    wage_frac_output = output_vs_wage.tot_wage.sum() / (va_2018_new.sum() + Z_2018_new.sum()).loc['221100']
    wage_to_va001 = wage_frac_output / va001_frac_output
    # impute_v001 has total wages rescaled to employee compensation using the same factor as
    # we empirically have for the electricity sector as a whole (employee comp ~= 2 x total wage)
    impute_v001 = output_vs_wage.frac_wage_of_output / wage_to_va001

    #impute_v001.plot.bar()

    # we impute V001 with our new estimates
    A_lit.loc['V001'] = impute_v001

    # 'Other electric power gen' spent > 100% on wages. This is unrealistic.
    # We set it to 30%. We don't use 'other' in our analysis, so this should have little impact
    A_lit.loc['V001', 'Other electric power generation'] = 0.30

    #1-va_frac

    # we readjust intermediate spending so all adds up to unity again
    # we see that solar reaches negative values
    interm_spending_frac_output = 1-va_frac-(A_lit.loc['V001'] - va001_frac_output)

    # we impute 2% for solar, while actually it is -4%
    impute_solar = 0.02
    A_lit.loc['V001', 'Solar electric power generation'] -= impute_solar - interm_spending_frac_output.loc['Solar electric power generation']

    # rebalance after imputation
    rebal_v001 = Z_va_2018_new.loc['V001', '221100'] / (A_lit.loc['V001'] * bea_electricity_output[str(io_year)]).sum()
    Z_v001 = (A_lit.loc['V001'] * bea_electricity_output[str(io_year)])  * rebal_v001
    A_lit.loc['V001'] = Z_v001 / bea_electricity_output[str(io_year)]

    interm_spending_frac_output = 1-va_frac-(A_lit.loc['V001'] - va001_frac_output)




    # readjust intermediate spending
    A_lit = pd.concat([(A_lit.drop(['V001', 'V002','V003']) * (interm_spending_frac_output / (1-va_frac))),
               A_lit.loc[['V001', 'V002','V003'], :]])

    # HERE DO SENSITIVITY ANALYSIS settings.SEED is set so the scramble can be recreated
    if settings.SENS_ANALYSIS_OPEX_COST_VECTORS == True:
        #To do min run
        #settings.SEED = 646
        #To do max run
        #settings.SEED = 584
        A_lit = help_sens.scramble_IO(A_lit,settings)

    assert_array_almost_equal(A_lit.sum().to_numpy(), np.ones(A_lit.shape[1], dtype=float))

    ### 4. We will impute the literature based values into the 2018 IO table. A number of conditions must hold:
    # i. Sum of inputs for disaggregated sectors equals inputs of original
    # ii. Total sum of input from disaggregated sectors equals input from original
    # iii. inputs to and from disaggregated dsectors equals original self-link
    # iV. All items > 0 (including final demand and value added)


    #### i. sum of inputs equals original inputs 

    # Scale Literature items so that output per electricity technology equals BEA numbers (recipes stay the same)
    Z_lit = A_lit.loc[:, bea_electricity_output.index].mul(bea_electricity_output[str(io_year)], axis=1)

    assert(Z_lit.sum().sum() - Z_va_2018_new.sum().loc['221100'] < 0.00000001)



    def ipfp_repeat(b, u, v, X, rowsum=False, n=2):
        a_n, b_n, m_n = ipfp(b, u, v, X, rowsum)
        for i in range(n):
            a_n, b_n, m_n = ipfp(b_n, u, v, X, rowsum)
        return a_n, b_n, m_n

    # the high number of zero-valued entries can lead to non-convergence. 
    # We set all rows where all but the default vectors are zero to the default vectors
    Z_ipfp = Z_lit.copy()
    lits = ['Hydroelectric power generation', 'Coal electric power generation', 'Gas electric power generation', \
                   'Solar electric power generation', 'Wind electric power generation', \
                   'Geothermal electric power generation', 'Biomass electric power generation']

    for ind in Z_ipfp.index:
        imp = True
        for tech in lits:
            if Z_ipfp.loc[ind, tech] > 0:
                imp = False
                break
        if imp:
            for tech in lits:
                Z_ipfp.loc[ind, tech] = bea_electricity_output.loc[tech, str(io_year)] * \
                                        A_lit.loc[ind, 'Other electric power generation'] / \
                                        (1 - va_frac + va001_frac_output - A_lit.loc['V001', 'Other electric power generation']) * \
                                        (1 - va_frac + va001_frac_output - A_lit.loc['V001', tech])       



    # Adjust matrix such that rows sum to BEA input figures, and columns sum to BEA output figures
    # additionally, we also set all zero-valued entries to a small value, so the ipfp converges
    # This cell laves the value added intact as they were
    small_value = 0.001
    va = ['V001', 'V002', 'V003']

    u = Z_va_2018_new['221100'].drop(va, axis=0).to_numpy()
    v = (bea_electricity_output[str(io_year)] * interm_spending_frac_output.loc[bea_electricity_output.index]).to_numpy()
    X = Z_ipfp.drop(va, axis=0).to_numpy()
    X[X == 0] = small_value
    b = np.ones(X.shape[1])

    a_n, b_n, m_n = ipfp_repeat(b, u, v, X, n=150)



    try:
        np.testing.assert_allclose(m_n.sum(axis=0), v)
    except AssertionError:
        return 1

    # adjustments can be made so that totals are within 99.8% of targets (relative difference < 0.002)
    try:
        np.testing.assert_allclose(m_n.sum(axis=1), u, rtol=7e-03 if imp_reduce < 1 else 1e-6)
    except AssertionError:
        return 1

    Z_lit_cols = pd.DataFrame(m_n, index=Z_lit.index.drop(va), columns=Z_lit.columns)
    # add value added again (unmuted)
    Z_lit_cols = pd.concat([Z_lit_cols, Z_lit.loc[va, :]])



    #### ii. sum of outputs = original output

    # The output of each sector is electricity
    # we simply split the output vectors to other sectors according to total output
    Z_lit_rows = pd.DataFrame(columns=Z_fd_2018_new.columns, index=Z_lit_cols.columns).fillna(0)
    for tech in Z_lit_rows.index:
        Z_lit_rows.loc[tech, :] = Z_fd_2018_new.loc['221100'] * \
                                (bea_electricity_output.loc[tech] / bea_electricity_output.sum())[str(io_year)]

    assert((Z_lit_rows.sum() - Z_fd_2018_new.loc['221100'] == 0).all)



    #### iii. Self links consistent

    Z_lit_core = pd.DataFrame(columns=Z_lit_cols.columns, index=Z_lit_cols.columns).fillna(1)

    # Literature self-links in BEA figures
    #Z_lit_cols.loc['221100', :]

    #Z_lit_rows.loc[:, '221100']

    # Adjust matrix such that rows sum to BEA input figures, and columns sum to BEA output figures

    u = Z_lit_rows.loc[:, '221100'].to_numpy()
    v = Z_lit_cols.loc['221100', :].to_numpy()
    X = Z_lit_core.to_numpy()
    b = np.ones(X.shape[1])

    a_n, b_n, m_n = ipfp_repeat(b, u, v, X, n=50)

    np.testing.assert_allclose(m_n.sum(axis=0), v)

    # accurate within 1% of row and column sums (relative difference  < 0.01)
    np.testing.assert_allclose(m_n.sum(axis=1), u, rtol=3e-3 if imp_reduce < 0.8 else 1e-7)

    Z_lit_core = pd.DataFrame(m_n, index=Z_lit_core.index, columns=Z_lit_core.columns)



    #### iv. All positive



    ### 5. Combine everything

    # build up the new Z matrix by inserting the new rows and columns in the right place
    # | old matrix        new column          old matrix |
    # | old matrix        new column          old matrix |
    # | new row           core                new row    |
    # | new row           core                new row    |
    # | old matrix        new column          old matrix |
    # | old matrix        new column          old matrix |
    # | old matrix        new column          old matrix |

    # Then add the value added and final demand submatrices in a similar fashion


    Z_2018_new_new = pd.concat([
        pd.concat([Z_2018_new.iloc[:5, :5], Z_lit_cols.iloc[:5, :], Z_2018_new.iloc[:5, 6:]], axis=1), 
        pd.concat([Z_lit_rows.iloc[:, :5], Z_lit_core, Z_lit_rows.iloc[:, 6:73]], axis=1), 
        pd.concat([Z_2018_new.iloc[6:, :5], Z_lit_cols.iloc[6:74, :], Z_2018_new.iloc[6:, 6:]], axis=1)],
        axis=0)

    fd_2018_new_new = pd.concat([fd_2018_new.iloc[:5, :], Z_lit_rows.iloc[:, 73:], fd_2018_new.iloc[6:, :]], axis=0)

    va_2018_new_new = pd.concat([va_2018_new.iloc[:, :5], Z_lit_cols.iloc[74:, :], va_2018_new.iloc[:, 6:]], axis=1)

    Z_2018_new_new.index = Z_2018_new_new.index.astype(str)
    va_2018_new_new.index = va_2018_new_new.index.astype(str)
    fd_2018_new_new.index = fd_2018_new_new.index.astype(str)

    Z_2018_new_new.columns = Z_2018_new_new.columns.astype(str)
    va_2018_new_new.columns = va_2018_new_new.columns.astype(str)
    fd_2018_new_new.columns = fd_2018_new_new.columns.astype(str)

    Z_fd_2018_new_new = pd.concat([Z_2018_new_new, fd_2018_new_new], axis=1)
    Z_va_2018_new_new = pd.concat([Z_2018_new_new, va_2018_new_new])
    Z_va_fd_2018_new_new = pd.concat([Z_fd_2018_new_new, va_2018_new_new])

    A_2018_new_new = Z_2018_new_new / Z_va_2018_new_new.sum()
    A_2018_new_new.loc[:, 'imports'] = 0

    # manual intervention biomass
    # we put 10% of biomass in agriculture: 5% in 111CA and 5% in 113FF (* 1-(value added frac))
    biomass_111CA = 0.05
    biomass_113FF = 0.05

    new = A_2018_new_new.loc[:, 'Biomass electric power generation']
    interm_frac = new.sum()

    new.loc['111CA'] = biomass_111CA
    new.loc['113FF'] = biomass_113FF

    interm_frac_after = new.sum()
    new *= (interm_frac / interm_frac_after)
    # sm is total output
    sm = Z_va_2018_new_new.loc[:, 'Biomass electric power generation'].sum()
    Z_2018_new_new.loc[:, 'Biomass electric power generation'] = new * sm

    Z_fd_2018_new_new = pd.concat([Z_2018_new_new, fd_2018_new_new], axis=1)
    Z_va_2018_new_new = pd.concat([Z_2018_new_new, va_2018_new_new])
    Z_va_fd_2018_new_new = pd.concat([Z_fd_2018_new_new, va_2018_new_new])

    A_2018_new_new = Z_2018_new_new / Z_va_2018_new_new.sum()
    A_2018_new_new.loc[:, 'imports'] = 0

    if imp_reduce == 1:
        Z_va_fd_2018_new_new.to_csv(settings.FL_data_sens + str(MODELRUN_ID) + '_Z_' + str(io_year) + '_inc_elec_split_lit.csv')

        #helper function
        Z_lit.to_csv(settings.FL_data_sens + str(MODELRUN_ID) + '_Z_lit.csv')
    else:
        Z_va_fd_2018_new_new.to_csv(settings.FL_data_sens + str(MODELRUN_ID) + '_Z_' + str(io_year) + '_inc_elec_split_lit_' + str(imp_reduce) + '.csv')

        #helper function
        Z_lit.to_csv(settings.FL_data_sens + str(MODELRUN_ID) + '_Z_lit_' + str(imp_reduce) + '.csv')
    

    L = np.linalg.inv(np.identity(A_2018_new_new.shape[0]) - A_2018_new_new)
    L_2018_new_new = pd.DataFrame(L, index=A_2018_new_new.index, columns=A_2018_new_new.columns)

    if imp_reduce == 1:
        L_2018_new_new.drop('imports', axis=1).to_csv(settings.FL_data_sens + str(MODELRUN_ID) + '_L_' + str(io_year) + '_inc_elec_split_lit.csv')
        A_2018_new_new.drop('imports', axis=1).to_csv(settings.FL_data_sens + str(MODELRUN_ID) + '_A_' + str(io_year) + '_inc_elec_split_lit.csv')
    else:
        L_2018_new_new.drop('imports', axis=1).to_csv(settings.FL_data_sens + str(MODELRUN_ID) + '_L_' + str(io_year) + '_inc_elec_split_lit_' + str(imp_reduce) + '.csv')
        A_2018_new_new.drop('imports', axis=1).to_csv(settings.FL_data_sens + str(MODELRUN_ID) + '_A_' + str(io_year) + '_inc_elec_split_lit_' + str(imp_reduce) + '.csv')
import numpy as np
import pandas as pd
import help_sens


def run(MODELRUN_ID, settings):

#get the years smoothing according to the sensitivity analysis that is being run
    # if not smoothing than smoothing years = 1
    if settings.SENS_ANALYSIS_SMOOTHING == True:
        years_smoothing = settings.YEARS_SMOOTHING_SENS
        smoothing = True
    else:
        years_smoothing = 1
        smoothing = False

    # get T and D values based on whether a sensitivity analysis is being run
    line_cost,thrice_line_cost_factor, thrice_line_opex_factor = help_sens.TandD_sens(settings)
    
    # Find cost in US of MW-mile of transmission

    #MW-miles
    NREL_trans = pd.read_csv(settings.FL_energy + 'transmission_capacity.csv', index_col=0)

    NREL_trans = NREL_trans.set_index(['years', 'scenario'])['capacity'].unstack()

    # We split this notebook in A: Capex, and B: Opex


    # A: Capex

    # Two reports by NREL
    # 1. Wind curtailment study (2017)	https://www.nrel.gov/docs/fy17osti/67240.pdf
    #    a. Average of Table 4: 1,384 per MW-mile
    #    b. Footnote 18: Other transmission cost estimates, 
    #        including transmission cost assumptions ranging from 
    #        about \\$900/MW-mile to over \\$3500/MW-mile used in the Wind Vision analysis,
    # 2. North america renewable integration (2021)	https://www.nrel.gov/docs/fy21osti/79224.pdf
    #    a. Data sources p. 14: Distance per mile costs were 
    #       estimated from the Eastern Interconnection Planning 
    #       Collaborative Phase 2 Report (EIPC 2015). Transmission costs 
    #       were \\$1,347\\$2,333 per MW-mile depending on the line voltage. 
    #       Spur lines connecting utility-scale solar and wind to the power 
    #       system cost \\$3,667/MWmile.

    # We pick 1384 / MW-mi

    # The 2017 cost index for power/comm structures (T50304A BEA) is 107.0; for 2018:110.8, 
    # so we increase the cost to 2018 USD:

    # 1384 * 110.8/107 = 1433

    # We follow Way et al. 2022 in their modelling approach. 
    # Instead of building new lines, the US will mostly need 
    # upgraded lines to deal with the high penetration of renewables, 
    # which can happen instead of planned maintenance and replacement 
    # (which is already in the system). Installing lines with 3x the 
    # capacity can be done at 1.37 the cost (Way et al SI p. 47).

    # 1. Add in between years
    # 2. Find change per year: this is new lines
    # 3. Divide by 2: This is the amount of miles of old line we are going to upgrade
    # 4. Multiply by 1384 x 1.37: this is the cost of upgrading the transmission network
    # 5. Multiply by (100/31), as transmission is 31% of total grid investment (distribution network is 69%)
    # 6. Smoothing 5 year
    # 7. Find the change in change per year: This is change in how many new lines are installed

    #1
    for yr in range(2011, 2050, 2):
        NREL_trans.loc[yr, :] = 0.5 * (NREL_trans.loc[yr-1, :] + NREL_trans.loc[yr+1, :])

    NREL_trans.sort_index(inplace=True)

    #2
    New_line_cost = NREL_trans.diff()
    New_line_cost.dropna(inplace=True)

    #3 #4
    New_line_cost = New_line_cost / 2 * thrice_line_cost_factor * line_cost

     # 5
    New_line_cost *= (100/31)

    #6
    # 5 year smoothing
    New_line_cost_sm = New_line_cost.copy()

    if smoothing == True:
        New_line_cost_sm = New_line_cost_sm.sort_index().reset_index()\
                        .rolling(years_smoothing).mean()
        New_line_cost_sm.rename(columns={'years': 'year'}, inplace=True)
        New_line_cost_sm = New_line_cost_sm[~New_line_cost_sm.year.isna()]
        New_line_cost_sm['year'] = New_line_cost_sm['year'].astype(int)
        New_line_cost_sm.set_index('year', inplace=True)

    New_line_cost_sm.to_csv(settings.FL_data_sens + str(MODELRUN_ID) + '_TandD_annual_invest_smooth' + str(smoothing) + '_18nov.csv')
    # total cost ~500 billion for 95% by 2035 scenario
    #New_line_cost_sm.sum()
    
    delta_line_invest = New_line_cost_sm.diff()
    delta_line_invest.dropna(inplace=True)
    #delta_line_invest.plot()

    delta_line_invest = delta_line_invest.loc[2020:, :]

    delta_line_invest.to_csv(settings.FL_data_sens + str(MODELRUN_ID) + '_TandD_annual_invest_delta_smooth' + str(smoothing) + '_18nov.csv')



    ## B Opex: what is the linear factor that employees have grown over the last 
    #          couple of years compared to the scenario change per MW-mile

    # network expansion from 2016 till 2018, and employment increase in same period
    #display(NREL_trans.loc[[2016, 2018], 'Reference'])

    diff = NREL_trans.loc[[2016, 2018], 'Reference'].values
    diff = diff[1] - diff[0]

    #print(diff)

    util_workers = pd.read_csv(settings.FL_preprocessed_data + 'util_workers_through_time.csv', index_col=0)

    gen = [221111, 221112, 221113, 221114, 221115, 221116, 221117, 221118]

    dist_trans = util_workers[util_workers.naics == 221100].groupby('yr').sum()['tot_emp'] -\
            util_workers[util_workers.naics.isin(gen)].groupby('yr').sum()['tot_emp']

    #print((dist_trans[18] - dist_trans[16]) / diff)

    #print(dist_trans[18] / 1.465204e+08)


    # delta employment total wobbles more than there is signal, so this won't be reliable

    #display(dist_trans)



    # instead, we use the same trick as in capex, and assume
    # employment growth x1.37 for every 3x grid expansion

    NREL_trans = NREL_trans.loc[2018:, :]

    # Additional workers are required proportional to the number of miles of line upgraded
    opex_factor = ( 
                    1.00 * (NREL_trans.loc[2018, :] - 0.5*(NREL_trans - NREL_trans.loc[2018, :])) + \
                    thrice_line_opex_factor * 0.5*(NREL_trans - NREL_trans.loc[2018, :])
                  ) / NREL_trans.loc[2018, :]

    # We use the final demand trick, where any change in output is added to final demand
    trans_dist = ['Electric bulk power transmission and control', 'Electric power distribution']
    Z_2018 = pd.read_csv(
        settings.FL_data_sens + str(MODELRUN_ID) + '_Z_2018_inc_elec_split_lit.csv',
        index_col=0,
    )

    #intermediate demand
    int_dem = Z_2018.loc[trans_dist].iloc[:, :-20].sum().sum()
    #print(int_dem)

    #final demand
    fin_dem = Z_2018.loc[trans_dist].iloc[:, -20:].sum().sum()
    #print(fin_dem)

    # total output
    tot_out = Z_2018.loc[:, trans_dist].sum().sum()
    #print(tot_out)

    frac_final_demand_of_total = fin_dem / (fin_dem + int_dem)

    opex_factor = ((opex_factor - 1) / frac_final_demand_of_total) + 1

    opex_factor.to_csv(settings.FL_data_sens + str(MODELRUN_ID) + '_TandD_opex_empl_factor_smooth' + str(smoothing) + '_18nov.csv')

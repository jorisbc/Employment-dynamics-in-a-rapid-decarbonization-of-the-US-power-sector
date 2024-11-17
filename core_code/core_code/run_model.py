import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import help_sens


def run(MODELRUN_ID, settings):
    # 0 SET PARAMETERS
    # 1 READ DATA
    # 2 DATA PREPROCESSING
    # 3 PLOTTING
    # 4 RUN MODEL PER SCENARIO PER YEAR

    # 0 SET PARAMETERS

    #get the years smoothing according to the sensitivity analysis that is being run
    # if not smoothing than smoothing years = 1
    if settings.SENS_ANALYSIS_SMOOTHING == True:
        years_smoothing = settings.YEARS_SMOOTHING_SENS
        smoothing = True
    else:
        years_smoothing = 1
        smoothing = False
            
    cost_name_prefix = ''
    if settings.SENS_ANALYSIS_COST:
        cost_name_prefix = help_sens.cost_sens(settings)
    
    io_year = settings.IO_YEAR_DEFAULT
    if settings.SENS_ANALYSIS_IO_YEAR:
        io_year = settings.IO_YEAR_SENS

    scens = ["95% by 2050", "95% by 2035", "Reference"]
    yrs = range(2020, 2050 - int((years_smoothing - 1)/2))
    capex_historic_baseline_years = [2019, 2020]
    if smoothing == False:
        yrs = range(2020, 2050)
        settings.FL_data_out += "_no_smoothing/"
    else:
        settings.FL_data_out += "/"
        if years_smoothing == 5:
            yrs = range(2021, 2050 - int((years_smoothing - 1)/2))


    # 1 READ DATA
    # A matrix 2018 domestic
    # if ((settings.SENS_ANALYSIS_OPEX_COST_VECTORS == False) and (settings.SENS_ANALYSIS_IO_YEAR == False)):
    #     A_sum_2018_dom = pd.read_csv(
    #         settings.FL_analysis
    #         + "Energy sector disaggregation/Scripts/output_lit/A_2018_inc_elec_split_lit.csv",
    #         index_col=0,
    #     )
    # else: # sensitivity analysis == True (either scrambled or other year, or both)
    A_sum_2018_dom = pd.read_csv(
        settings.FL_data_sens + str(MODELRUN_ID) + '_'
        + "A_" + str(io_year) + "_inc_elec_split_lit.csv",
        index_col=0,
    )

    # Z matrix 2018 domestic for total output numbers and final demand
    Z_2018 = pd.read_csv(
        settings.FL_data_sens + str(MODELRUN_ID) + '_'
        + "Z_" + str(io_year) + "_inc_elec_split_lit.csv",
        index_col=0,
    )


    # capex vectors
    investment_weights = pd.read_excel(
        settings.FL_data + "Synthetic IO sectors energy technologies alt.xlsx",
        sheet_name="Cost vectors with fuel cost",
        index_col=[0, 1],
        header=[0, 1],
        skipfooter=2,
    )
    
    # employment per industry output
    prefix_emp_per_mill = settings.OCC_DEFAULT_PREFIX
    if settings.SENS_ANALYSIS_OCC:
        prefix_emp_per_mill = settings.OCC_SENS_PREFIX
    emp_per_million_ouput_sum_2018 = pd.read_csv(
        settings.FL_data_sens + prefix_emp_per_mill + "emp_per_million_output_sum_2018.csv", index_col=0
    )

    # occupation names
    occ_names = pd.read_csv(
        settings.FL_preprocessed_data + "occ_names_bls_minor_major.csv", index_col=0
    )

    # percentage import per industry (we only use domestic jobs)
    import_perc = pd.read_excel(settings.FL_data + "Data_out_IO/2018-import_perc_sum.xlsx", index_col=0)


    # settings of colors for figures
    colornames = [
    "black","yellow","blue","magenta","limegreen","steelblue","pink","teal","tomato","darkblue",
    "crimson","aqua","darkviolet","orange","olive","brown","grey","lime","fuchsia","peachpuff",
    "g","y",]

    colorscheme = pd.DataFrame(
        occ_names["class_major"].unique(), columns=["classification"]
    )
    colorscheme["color"] = colornames
    occ_names = (
        occ_names.reset_index()
        .merge(colorscheme, right_on="classification", left_on="class_major")
        .set_index("OCC_CODE")
    )

    
    # setup output dataframes
    ## emp_invest
    emp_invest = occ_names.loc[emp_per_million_ouput_sum_2018.index].copy()
    emp_invest.drop("index", axis=1, inplace=True)

    ## new_jobs
    ## copy emp_invest index and occupation titles and delete all other columns
    new_jobs = occ_names.loc[emp_invest.index].copy()
    new_jobs.drop("index", axis=1, inplace=True)
    new_jobs.drop(new_jobs.iloc[:,1:15], axis=1, inplace=True)

    # rename sectors
    title_capitals = {
        "Biomass electric power generation": "Biomass Electric Power Generation",
        "Electric bulk power transmission and control": "Electric Bulk Power Transmission and Control",
        "Electric power distribution": "Electric Power Distribution",
        "Geothermal electric power generation": "Geothermal Electric Power Generation",
        "Hydroelectric power generation": "Hydroelectric Power Generation",
        "221200": "Natural Gas Distribution",
        "Nuclear electric power generation": "Nuclear Electric Power Generation",
        "Other electric power generation": "Other Electric Power Generation",
        "Solar electric power generation": "Solar Electric Power Generation",
        "221300": "Water, Sewage and Other Systems",
        "Wind electric power generation": "Wind Electric Power Generation",
        "Coal electric power generation": "Coal Electric Power Generation",
        "Gas electric power generation": "Gas Electric Power Generation",
    }
    title_lower = dict((v, k) for k, v in title_capitals.items())

    emp_per_million_ouput_sum_2018.rename(title_lower, axis=1, inplace=True)
    emp_per_million_ouput_sum_2018.drop("GFGD", axis=1, inplace=True)


    # read scenario data from folder including sensitivity analysis data
    folder = settings.FL_data_sens + "cost_trajectories/"
    scen_data = pd.read_csv(folder + cost_name_prefix + "alldata_io_tech_nosmoothing.csv", index_col=0)
    scen_data.rename(columns={"years": "year"}, inplace=True)

    tech_names = {
        "Bio": "Biomass electric power generation",
        "Coal": "Coal electric power generation",
        "Solar": "Solar electric power generation",
        "Gas": "Gas electric power generation",
        "Geo": "Geothermal electric power generation",
        "Hydro": "Hydroelectric power generation",
        "Nuclear": "Nuclear electric power generation",
        "Wind": "Wind electric power generation",
    }
    oth = ["Other electric power generation", "GFGD"]


    # ### 2 PREPROCESSING STEPS
    # 0. deflate to 2018-USD
    # 1. five year smoothing
    # 2. add baseline cost
    # 3. find capex 'shock', i.e. difference from year before
    # 4. separate distributed solar (rooftop, commercial) from utility solar
    # 5. calculate factor from USD cost to USD output in 2018
    # 6. add T&D cost (capex and opex)

    # 2.0 delfate to 2018 USD: Everything is in 2019 USD except coal/gas fuel which are in 2018 USD already:
    # 0: inflation 2018-2019 was 1.8%
    # 1: not coal/gas: investment, annual.capex, annual.FOM, annual.VOM, annual.Fuel, annual.opex: divide by 1.018
    # 2: coal/gas: a) investment, annual.capex, annual.VOM, annual.FOM: divide 1.018 for gas/coal,
    #              b) annual.opex = sum(annual.VOM, .FOM, .Fuel)

    # 0
    scen_data.loc[
        ~scen_data.technology.isin(["Gas", "Coal"]),
        [
            "investment",
            "annual.capex",
            "annual.FOM",
            "annual.VOM",
            "annual.Fuel",
            "annual.opex",
        ],
    ] /= 1.018

    # 1
    scen_data.loc[
        scen_data.technology.isin(["Gas", "Coal"]),
        ["investment", "annual.capex", "annual.FOM", "annual.VOM"],
    ] /= 1.018
    # 2
    scen_data.loc[
        scen_data.technology.isin(["Gas", "Coal"]), "annual.opex"
    ] = scen_data.loc[
        scen_data.technology.isin(["Gas", "Coal"]),
        ["annual.FOM", "annual.VOM", "annual.Fuel"],
    ].sum(
        axis=1
    )

    scen_data.loc[
        scen_data.technology.isin(["Coal", "Gas"]) & (scen_data.year < 2019), "annual.opex"
    ] = np.nan

    # 2.1 5 year smoothing
    scen_sm = scen_data.copy()
    if smoothing == True:
        scen_sm = (
            scen_sm.sort_values("year")
            .groupby(["technology", "scenario"])
            .rolling(years_smoothing)
            .mean()
            .reset_index()
            .drop("level_2", axis=1)
        )
        scen_sm = scen_sm[~scen_sm.year.isna()]
        scen_sm["year"] = scen_sm["year"].astype(int)

    # 2.2 Add baseline cost
    baseline_capex = (
        scen_data[scen_data.year.isin(capex_historic_baseline_years)]
        .groupby(["technology", "scenario"])["annual.capex"]
        .mean()
    )


    scen_sm["delta_capex"] = scen_sm["annual.capex"]

    # set baseline to 2020
    for tech in set(scen_sm.technology):
        for scen in set(scen_sm.scenario):
            scen_sm.loc[
                (scen_sm.technology == tech)
                & (scen_sm.scenario == scen)
                & (scen_sm.year == 2020),
                "delta_capex",
            ] = baseline_capex.loc[(tech, scen)]

    # 2.3 Calculate capex shock
    # capex shock is difference with previous year capex spending
    capex_shock = scen_sm[["year", "technology", "scenario", "delta_capex"]]
    capex_shock = capex_shock[capex_shock.year > 2019]
    #capex_shock["delta_capex"] = capex_shock.groupby(["technology", "scenario"])[
    #    "delta_capex"
    #].diff()
    #capex_shock = capex_shock[~(capex_shock.year == 2020)]

    if smoothing == False:
        capex_shock = (
            capex_shock.sort_values("year")
            .groupby(["technology", "scenario"])
            .rolling(2)
            .mean()
            .reset_index()
            .drop("level_2", axis=1)
        )
        capex_shock = capex_shock[~capex_shock.year.isna()]
        capex_shock["year"] = capex_shock["year"].astype(int)
        capex_shock = capex_shock[["year", "technology", "scenario", "delta_capex"]]

    # 2.4 Separate solar rooftop from solar utility
    folder = settings.FL_data_sens + "cost_trajectories/"
    all_data = pd.read_csv(folder + cost_name_prefix + "alldata_scenario_tech.csv", index_col=0)

    year_set_fact = [2020]
    calc_fact = all_data[all_data.years.isin(year_set_fact)]

    # solar distributed (rooftop + commercial firms with own solar) is treated separately
    tech_io = {
        "Biomass": "Biomass electric power generation",
        "CSP": "Solar electric power generation",
        "Coal": "Coal electric power generation",
        "Gas CC": "Gas electric power generation",
        "Gas CC CCS": "Gas electric power generation",
        "Gas CT": "Gas electric power generation",
        "Geothermal": "Geothermal electric power generation",
        "Hydro": "Hydroelectric power generation",
        "Nuclear": "Nuclear electric power generation",
        "PHS": "Hydroelectric power generation",
        "PV dist": "Solar dist",
        "PV utility": "Solar electric power generation",
        "Wind off": "Wind electric power generation",
        "Wind on": "Wind electric power generation",
    }
    calc_fact.technology.replace(tech_io, inplace=True)


    # 2.5 calculate cost to output factor
    calc_fact = calc_fact.groupby(["years", "technology", "scenario"]).sum()


    # deflate to 2018 USD: Everything is in 2019 USD except coal/gas fuel which are in 2018 USD already:
    # 0: inflation 2018-2019 was 1.8%
    # 1: not coal/gas: investment, annual.capex, annual.FOM, annual.VOM, annual.Fuel, annual.opex: divide by 1.018
    # 2: coal/gas: a) investment, (annual).capex, (annual.)VOM, (annual.)FOM: divide 1.018 for gas/coal,
    #              b) annual.opex = sum(annual.VOM, .FOM, .Fuel)

    # 0
    calc_fact.loc[
        (
            slice(None),
            ~calc_fact.index.isin(
                ["Gas electric power generation", "Coal electric power generation"],
                level="technology",
            ),
        ),
        [
            "investment",
            "CAPEX",
            "FOM",
            "VOM",
            "Fuel",
            "annual.capex",
            "annual.FOM",
            "annual.VOM",
            "annual.Fuel",
            "annual.opex",
        ],
    ] /= 1.018
    # 1
    calc_fact.loc[
        (
            slice(None),
            ["Gas electric power generation", "Coal electric power generation"],
            slice(None),
        ),
        [
            "investment",
            "CAPEX",
            "FOM",
            "VOM",
            "annual.capex",  #'Fuel', 'annual.Fuel', 'annual.opex',
            "annual.FOM",
            "annual.VOM",
        ],
    ] /= 1.018
    # 2
    calc_fact.loc[
        (
            slice(None),
            ["Gas electric power generation", "Coal electric power generation"],
            slice(None),
        ),
        "annual.opex",
    ] = calc_fact.loc[
        (
            slice(None),
            ["Gas electric power generation", "Coal electric power generation"],
            slice(None),
        ),
        ["annual.FOM", "annual.VOM", "annual.Fuel"],
    ].sum(
        axis=1
    )


    elec_output = pd.read_excel(
        settings.FL_data + "BEA/GrossOutput_2021.xlsx",
        sheet_name="UGO305-A",
        usecols="B,AA",
        index_col=0,
        skiprows=29,
        nrows=8,
        header=None,
    )

    elec_output.columns = year_set_fact
    calc = calc_fact.loc[(year_set_fact, slice(None), "95% by 2035"), "annual.opex"]
    calc.index = calc.index.droplevel([0, 2])

    # factor calc:
    factor = (calc / 1000) / elec_output[2020]
    factor = factor.drop(
        [
            "Batteries",
            "Fossil fuel electric power generation",
            "Other electric power generation",
        ]
    )

    factor["Coal electric power generation"] = (
        calc.loc[["Coal electric power generation", "Gas electric power generation"]].sum()
        / 1000
    ) / elec_output.loc["Fossil fuel electric power generation", year_set_fact]
    factor["Gas electric power generation"] = (
        calc.loc[["Coal electric power generation", "Gas electric power generation"]].sum()
        / 1000
    ) / elec_output.loc["Fossil fuel electric power generation", year_set_fact]

    fact = pd.DataFrame(factor, columns=year_set_fact)
    fact = fact[2020].drop("Solar dist")


    # 2.6 Add T&D cost
    # add T&D capex spending (see other file)
    trans_dist_capex = pd.read_csv(
        settings.FL_data_sens + str(MODELRUN_ID) + "_TandD_annual_invest_smooth" + str(smoothing) + "_18nov.csv"
    )

    trans_dist_capex = trans_dist_capex[trans_dist_capex.year >= 2020]
    trans_dist_capex["technology"] = "Trans&dist"
    trans_dist_capex = trans_dist_capex.set_index(["year", "technology"]).stack() / 1000
    trans_dist_capex.index.set_names("scenario", level=-1, inplace=True)
    trans_dist_capex = trans_dist_capex.reset_index().rename(columns={0: "delta_capex"})

    # merge on capex_shock
    capex_shock = pd.concat([capex_shock, trans_dist_capex])
    capex_shock.sort_values(by=["technology", "scenario", "year"], inplace=True)
    capex_shock.reset_index(drop=True, inplace=True)

    # T&D opex spending factor
    trans_dist_opex_fact = pd.read_csv(
        settings.FL_data_sens + str(MODELRUN_ID) +  "_TandD_opex_empl_factor_smooth" + str(smoothing) + "_18nov.csv"
    )

    trans_dist_opex_fact = trans_dist_opex_fact.set_index("years").stack()
    trans_dist_opex_fact.index.set_names("scenario", level=-1, inplace=True)

    if settings.SENS_ZERO_DOUBLE_SOLAR_WIND:
        if settings.SOLAR_WIND_SENS_OPTION == 'zero':
            capex_shock.loc[capex_shock.technology == 'Wind', 'delta_capex'] = 0
            capex_shock.loc[capex_shock.technology == 'Solar', 'delta_capex'] = 0
            scen_sm.loc[scen_sm.technology == 'Wind', 'annual.opex'] = 0
            scen_sm.loc[scen_sm.technology == 'Solar', 'annual.opex'] = 0
        elif settings.SOLAR_WIND_SENS_OPTION == 'double':
            scen_sm.loc[scen_sm.technology == 'Wind', 'annual.opex'] *= 2
            scen_sm.loc[scen_sm.technology == 'Solar', 'annual.opex'] *= 2
        else:
            # raise error
            raise ValueError('SOLAR_WIND_SENS_OPTION value is not recognized: can be "zero" or "double"') 


    # 3 Plotting
    # two charts of electricity generation: in MWh and in opex cost
    gen_graph = scen_sm.loc[
        scen_sm.scenario == scens[1], ["technology", "year", "generation", "annual.opex"]
    ]
    gen_graph = gen_graph.set_index(["year", "technology"]).unstack()

    gen_graph["generation"].plot.line(figsize=(7.5, 5)).legend(
        loc="center left", bbox_to_anchor=(1.0, 0.5)
    )
    plt.tight_layout()
    
    # Only create figure if settings.CREATE_FIGURES == True:
    if settings.CREATE_FIGURES == True:
        plt.savefig(
            settings.FL_fig_folder
            + str(MODELRUN_ID)
            + "_generation_scenario_"
            + scens[1]
            + "_smooth_"
            + str(smoothing)
            + "_2020-2050.png"
        )

    gen_graph["annual.opex"].plot.line(figsize=(7, 5)).legend(
        loc="center left", bbox_to_anchor=(1.0, 0.5)
    )
    plt.tight_layout()

    # Only create figure if settings.CREATE_FIGURES == True:
    if settings.CREATE_FIGURES == True:    
        plt.savefig(
            settings.FL_fig_folder
            + str(MODELRUN_ID)
            + "_generation_cost_scenario_"
            + scens[1]
            + "_smooth_"
            + str(smoothing)
            + "_2020-2050.png"
        )

    # chart of change in capex spending per year
    for scen in scens:
        delta_capex = (
            capex_shock[capex_shock.scenario == scen]
            .drop("scenario", axis=1)
            .set_index(["year", "technology"])
            .unstack()
        )

        delta_capex["delta_capex"].plot.line(figsize=(7.5, 5)).legend(
            loc="center left", bbox_to_anchor=(1.0, 0.5)
        )
        plt.tight_layout()

        # Only create figure if settings.CREATE_FIGURES == True:
        if settings.CREATE_FIGURES == True:
            plt.savefig(
                settings.FL_fig_folder
                + str(MODELRUN_ID)
                + "_capex_per_year_scen_"
                + scen
                + "_smooth_"
                + str(smoothing)
                + "_2020-2050.png"
            )


    # 4 RUN FOR EACH YEAR: OPEX AND CAPEX

    # 4.1 Help functions
    def adjust_io_elec(A, s_gen, yr, scenario, yr_b=None):
        # gen gives electricity mix
        gen = s_gen.loc[(yr, scenario)]

        # if we use a base year, we also calculate the increase in electricity use
        elec_increase_fact = 1.0
        if yr_b is not None:
            gen_b = s_gen.loc[(yr_b, scenario)]
            gen_b.index = gen_b.index.map(tech_names)
            elec_increase_fact = gen.sum() / gen_b.sum()

        # if settings.SENS_IMPEXP == True: # update IO matrix if we are doing a sensitivity analysis on import/export
        #     A, _, _ = help_sens.imp_sens(MODELRUN_ID, yr, settings)

        # A_tot has the total fraction of spending per industry on electricity sectors
        A_tot = A.loc[list(tech_names.values()), :].sum()
        gen_frac = gen.copy()
        gen_frac /= gen_frac.sum()
        gen_frac *= elec_increase_fact
        # We take the outer product of the electricity mix and spending per industry
        A_new_elec = pd.DataFrame(
            np.outer(gen_frac, A_tot), index=gen_frac.index, columns=A_tot.index
        )

        A_new = A.copy()

        # adjust so that the increase in elec spending does not affect it to sum to unity
        #             a_elec_change = A_new_elec.sum(axis=0) * (elec_increase_fact-1.0) # change in fraction spent on elec
        #             other_a_sum = A_new.drop(A_new_elec.index, axis=0).sum(axis=0) # sum of all other spending fractions
        #             other_a_sum = other_a_sum / (other_a_sum + a_elec_change) # adjustment needed
        #             A_new *= other_a_sum # adjust A_new

        # We insert the new electricity generation sectors' inputs in the A matrix
        A_new.loc[A_new_elec.index, :] = A_new_elec
        A_new.loc[:, oth] = A.loc[:, oth]
        A_new.loc[
            "Other electric power generation"
        ] = 0  # All electricity is accounted for, there is no 'other'

        # We add an imports sector to A in order to calculate the L matrix
        A_new_imp = A_new.copy()
        A_new_imp["imports"] = 0
        L_new = np.linalg.inv(np.identity(A_new_imp.shape[0]) - A_new_imp)
        L_new = pd.DataFrame(L_new, index=A_new_imp.index, columns=A_new_imp.columns)

        return L_new.iloc[:, :-1], A_new


    def extract_fd_from_Z(X, L, fd, idx):
        # X = L fd
        C = L.loc[idx].drop(idx, axis=1) @ fd.drop(idx)
        tot = L @ fd
        fd_idx = np.linalg.solve(np.linalg.inv(L.loc[idx, idx]), X.loc[idx] - C)
        fd.loc[idx] = fd_idx
        return fd


    def fd_opex(fd, s_gen, yr, scenario, L):  # , yr_b = None):
        gen = s_gen.loc[(yr, scenario)]

        # gen_base = s_gen.loc[(2020, scenario)]

        # gen is new output in dollars, so part of the Z matrix
        # X = L fd. We need to find the elec values of fd that reproduce the new Z values
        fd_new = extract_fd_from_Z(
            gen.loc[tech_names.values()], L, fd.copy(), tech_names.values()
        )

        # no_change = fd_imi_2018.loc[tech_names.values()].sum() / delta_fd_elec.sum()
        # fd.loc[delta_fd_elec.index] = delta_fd_elec * no_change
        # return fd.copy()
        return fd_new


    def L_times_delta_fd(L_b, L_e, fd_b, fd_e):
        if fd_b.size == 1:
            x_b = L_b * fd_b
            x_e = L_e * fd_e
        else:
            x_b = L_b @ fd_b
            x_e = L_e @ fd_e
        return x_e, x_b, x_e - x_b


    def adjust_io_and_fd_inverse_calc(A, fd, s_gen, yr_b, yr_e, scenario):
        opex_output = pd.DataFrame(index=A.index, columns=tech_names.values())
        # final demand base year and end year
        L_b, A_b = adjust_io_elec(A, s_gen, yr_b, global_scen)
        L_e, A_e = adjust_io_elec(A, s_gen, yr_e, global_scen, yr_b=yr_b)

        fd_b = fd_opex(fd, s_gen, yr_b, scenario, L_b)
        fd_e = fd_opex(fd, s_gen, yr_e, scenario, L_e)

        # for tech in tech:
        # a, b can be used as sanity check to see if output difference matches electricity generation change
        a, b, opex_output.loc[:, "indirect"] = L_times_delta_fd(
            L_b.drop(tech_names.values(), axis=1),
            L_e.drop(tech_names.values(), axis=1),
            fd_b.drop(tech_names.values()),
            fd_e.drop(tech_names.values()),
        )

        for tech in tech_names.values():
            c, d, opex_output.loc[:, tech] = L_times_delta_fd(
                L_b[tech], L_e[tech], fd_b.loc[tech], fd_e.loc[tech]
            )
            a += c
            b += d

        return a, b, L_b, L_e, A_b, A_e, opex_output, fd_b, fd_e


    gen_2018 = scen_sm.loc[
        (scen_sm.year == 2018) & (scen_sm.scenario == "Reference"), "generation"
    ].sum()
    trans_dist = [
        "Electric bulk power transmission and control",
        "Electric power distribution",
    ]

        # accomodate the sensitivity analysis for CAPEX I/O matrix
    if settings.SENS_ANALYSIS_CAPEX_COST_VECTORS == True: # or settings.SENS_ANALYSIS_OPEX_COST_VECTORS == True:
        #get an array of sub-runs for sensitivity analysis involving sub_runs
        subrun_array = np.random.randint(low=0, high=1000, size = settings.SENS_ANALYSIS_SUBRUNS)

        #save the sub_runs seed numbers in the data_out folder
        df_subrun_array = pd.DataFrame(subrun_array)
        df_subrun_array.to_csv(settings.FL_data_out + str(MODELRUN_ID) + "_capex_subruns.csv", index=False)
    else:
        # if a seed exists use it (enables a replication of a model run that uses a random number generator)
        # put whatever seed value is used into an array so it can be iterated in loop below
        if settings.USE_SEED == True:
            subrun_array = np.full(1, settings.SEED)
        else:
            subrun_array = np.full(1, 0)

    # 4.2 Run for each scenario and every year
    for sub_run_num in subrun_array:
        #count the sub_runs for naming output files 
        # (for SENS_ANALYSIS_OPEX_COST_VECTORS the loop is in Energy_occ_io_sensitivity)
        settings.sub_run_id += 1
        # For sensitivity analysis
        if settings.SENS_ANALYSIS_CAPEX_COST_VECTORS == True or settings.USE_SEED == True:
            #use the random number generated for the sub_run as the seed
            # if USE_SEED == True then use the SEED to run the scramble_IO as the user is trying to replicate a subrun from a sensitivity analysis
            if settings.SENS_ANALYSIS_CAPEX_COST_VECTORS == True:
                settings.SEED = sub_run_num 
            
            #To do min run
            #settings.SEED = 128
            #To do max run
            #settings.SEED = 390
            investment_weights_modelrun = help_sens.scramble_IO(investment_weights,settings)
        else:
            investment_weights_modelrun = investment_weights.copy()

        for global_scen in scens:
            for yr_b in yrs:
                yr_e = yr_b + 1

                if settings.SENS_IMPEXP == True: # update IO matrix if we are doing a sensitivity analysis on import/export
                    import_perc = help_sens.imp_sens(MODELRUN_ID, yr_e, settings) #A_sum_2018_dom, Z_2018, 
                    import_perc_prev = help_sens.imp_sens(MODELRUN_ID, yr_b, settings) #A_sum_2018_dom, Z_2018, 
                    # A matrix 2018 domestic with reduced imports: 2030: 0.9x, 2035: 0.7x, 2040: 0.5x
                else:
                    import_perc = import_perc_prev = import_perc

                # 1. PREP GENERATION SCENARIO
                gen_scen = scen_sm[
                    ["year", "scenario", "technology", "generation", "annual.opex"]
                ].copy()

                gen_scen_batt = gen_scen[(gen_scen.technology == "Batteries")]
                # we add batteries back later: we do not have opex spending pattern, only capex
                # We assume batteries opex spending reflects capex spending
                gen_scen = gen_scen[~(gen_scen.technology == "Batteries")]
                gen_scen = gen_scen[gen_scen.year > 2019]

                gen_scen = gen_scen.set_index("technology")
                gen_scen = gen_scen.rename(index=tech_names)
                gen_scen["fact"] = fact

                gen_scen = (
                    gen_scen.reset_index().groupby(["year", "scenario", "technology"]).mean()
                )
                gen_scen["annual.opex"] *= 1000

                gen_scen[
                    "annual.opex"
                ] /= gen_scen.fact  # factor to transform cost into revenue

                # setup dataframes to store sensitivity analysis outputs
                #-	Maximum net new jobs since 2020
                #new_jobs = pd.DataFrame()#{'OCC_CODE':[],'OCC_TITLE':[],'NEW_JOBS':[],'MAX_NEW_JOBS':[],'MAX_NEW_JOBS_YEAR':[]})
                # Setting index for this Dataframe
                #set_index = new_jobs.set_index(['OCC_CODE'], inplace=True)
            
                # Display DataFrame
                #print("DataFrame 1:\n",new_jobs)
                # add solar distrubited (IO only has solar utility)
                for yr in [yr_b, yr_e]:
                    for col in ["generation", "annual.opex"]:
                        # dist_solar_factor = pv_dist_cost / pv_utility_cost
                        dist_solar_factor = (
                            all_data.loc[
                                (all_data.years == yr)
                                & (all_data.scenario == global_scen)
                                & (all_data.technology == "PV dist"),
                                col,
                            ].values[0]
                            / all_data.loc[
                                (all_data.years == yr)
                                & (all_data.scenario == global_scen)
                                & (all_data.technology == "PV utility"),
                                col,
                            ].values[0]
                        )

                        gen_scen.loc[
                            (yr, global_scen, "Solar electric power generation"), col
                        ] *= (1 + dist_solar_factor)

                # If we don't use the scaling factor, we could add a profit margin to the cost and use that instead
                # gen_scen['annual.opex'] *= 1.46 # add profit margin to cost (cost + profit = revenue)

                # 2. CALCULATE GENERATION SCENARIO
                # 2.1 We calculate (domestic) final demand from the domestic Z matrix
                fd_imi_2018 = Z_2018.iloc[:-4, -20:].sum(axis=1)

                (
                    a,
                    b,
                    L_b,
                    L_calc,
                    A_b,
                    A_calc,
                    opex_output,
                    fd_b,
                    fd_e,
                ) = adjust_io_and_fd_inverse_calc(
                    A_sum_2018_dom.copy(),
                    fd_imi_2018.copy(),
                    gen_scen["annual.opex"].copy() / 1000000,
                    yr_b,
                    yr_e,
                    global_scen,
                )

                opex_output_year = opex_output / (yr_e - yr_b) * 1000000

                # the indirect impact is originally 136% of total impact
                total_fact = (
                    Z_2018.loc[tech_names.values()].sum().sum()
                    / fd_imi_2018.loc[tech_names.values()].sum()
                )

                # the indirect impact now is a certain factor of total impact
                opex_output_total = opex_output_year.loc[
                    :, tech_names.values()
                ].sum(axis=1)
                zeros = (opex_output_total == 0)
                opex_output_total.loc[zeros] = 1
                indir_fact = opex_output_year.indirect / opex_output_total
                indir_fact.loc[zeros] = 0  # some division by 0 if total output change = 0
                # We assign the indirect impact to the different technologies
                opex_output_year = opex_output_year.multiply(indir_fact + 1, axis=0)

                opex_output_year.drop(["indirect"], axis=1, inplace=True)

                # Increasing electricity use leads to higher output for transmission and distribution sectors
                fd_e1 = fd_e.copy()
                fd_e2 = fd_e.copy()
                fd_e1.loc[trans_dist] *= trans_dist_opex_fact.loc[yr_b, global_scen]
                fd_e2.loc[trans_dist] *= trans_dist_opex_fact.loc[yr_e, global_scen]

                opex_output_year["Transmission and distribution"] = (
                    L_calc @ fd_e2 - L_calc @ fd_e1
                ) * 1000000

                # 2.2 capex investments
                # 2.2.1 prep capex weights
                investment_w = investment_weights_modelrun.copy()
                investment_w.index = investment_w.index.droplevel(0)
                investment_w.index = investment_w.index.astype(str)

                #### adjust investment weights with A matrix weights per electricity tech

                # fraction of utilities spending on direct gas (221200) and water and waste mngment (221300)
                frac_221200 = (
                    Z_2018.iloc[5:18].sum(axis=1).loc["221200"]
                    / Z_2018.iloc[5:18].sum(axis=1).sum()
                )
                frac_221300 = (
                    Z_2018.iloc[5:18].sum(axis=1).loc["221300"]
                    / Z_2018.iloc[5:18].sum(axis=1).sum()
                )

                # only take part of investment that is spent locally
                investment_w_curr = investment_w.multiply(
                    1.0 - import_perc["import_perc"].loc[investment_w.index], axis=0
                )
                investment_w_prev = investment_w.multiply(
                    1.0 - import_perc_prev["import_perc"].loc[investment_w.index], axis=0
                )

                # split up spending on utilities sector in electricity and gas and waste management
                investment_w_curr.loc["221200"] = frac_221200 * investment_w_curr.loc["22"]
                investment_w_curr.loc["221300"] = frac_221300 * investment_w_curr.loc["22"]
                investment_w_curr.loc["22"] *= 1 - frac_221200 - frac_221300
                investment_w_prev.loc["221200"] = frac_221200 * investment_w_prev.loc["22"]
                investment_w_prev.loc["221300"] = frac_221300 * investment_w_prev.loc["22"]
                investment_w_prev.loc["22"] *= 1 - frac_221200 - frac_221300

                # further split it into different electricity generation sectors
                elec_dist = (A_calc.iloc[5:16] / A_calc.iloc[5:16].sum()).mean(axis=1)
                inv_elec = pd.DataFrame(
                    np.outer(elec_dist.to_numpy(), investment_w_curr.loc["22"].to_numpy()),
                    index=elec_dist.index,
                    columns=investment_w_curr.loc["22"].index,
                )
                inv_elec = pd.DataFrame(
                    np.outer(elec_dist.to_numpy(), investment_w_prev.loc["22"].to_numpy()),
                    index=elec_dist.index,
                    columns=investment_w_prev.loc["22"].index,
                )

                investment_w_curr = pd.concat([investment_w_curr, inv_elec])
                investment_w_curr.drop("22", inplace=True)

                investment_w_prev = pd.concat([investment_w_prev, inv_elec])
                investment_w_prev.drop("22", inplace=True)

                # Pumped hydro investments are added to 'hydro', so we do not
                # need a separate explicit 'pumped hydro' sector
                inv = investment_w_curr.drop(["Pumped hydro"], axis=1)
                inv_prev = investment_w_prev.drop(["Pumped hydro"], axis=1)

                inv = inv.loc[:, pd.IndexSlice[:, "capex"]]
                inv_prev = inv_prev.loc[:, pd.IndexSlice[:, "capex"]]

                # 2.2.2: add capex scenario
                cap_year_shock = (
                    capex_shock.set_index(["year", "scenario", "technology"])
                    .stack()
                    .unstack("technology")
                )
                cap_2030_shock = cap_year_shock.loc[2030].groupby(level="scenario").sum()
                cap_year_shock_prev = cap_year_shock.loc[[yr for yr in range(yr_b, yr_e)]].groupby(
                    level='scenario').sum()
                cap_year_shock = cap_year_shock.loc[[yr + 1 for yr in range(yr_b, yr_e)]].groupby(
                    level='scenario').sum()


                # Note that coal and nuclear opex are left out
                shck_index = {
                    "Bio": "Biomass",
                    "Gas": "Natural gas",
                    "Batteries": "Battery storage",
                    "Geo": "Geothermal",
                    "Trans&dist": "T&D",
                }

                cap_year_shock.rename(columns=shck_index, inplace=True)
                cap_year_shock_prev.rename(columns=shck_index, inplace=True)
                cap_2030_shock.rename(columns=shck_index, inplace=True)

                for col in inv.columns.get_level_values(0):
                    inv.loc[:, (col, "capex")] *= cap_year_shock.loc[global_scen, col] * 1000
                for col in inv_prev.columns.get_level_values(0):
                    inv_prev.loc[:, (col, "capex")] *= cap_year_shock_prev.loc[global_scen, col] * 1000
                

                if settings.SENS_IMPEXP == True: # update IO matrix if we are doing a sensitivity analysis on import/export
                    inv_exp_wind, inv_exp_solar = help_sens.exp_sens(MODELRUN_ID, yr_e, cap_2030_shock, import_perc["import_perc"], investment_w.copy(), global_scen, settings)
                    # element-wise multiplication for Wind and Solar only
                    if inv_exp_wind is not None:
                        inv.loc[inv_exp_wind.index, ("Wind", "capex")] += inv_exp_wind
                        inv.loc[inv_exp_solar.index, ("Solar", "capex")] += inv_exp_solar

                    inv_exp_wind_prev, inv_exp_solar_prev = help_sens.exp_sens(MODELRUN_ID, yr_b, cap_2030_shock, import_perc_prev["import_perc"], investment_w, global_scen, settings)
                    # element-wise multiplication for Wind and Solar only
                    if inv_exp_wind_prev is not None:
                        inv_prev.loc[inv_exp_wind_prev.index, ("Wind", "capex")] += inv_exp_wind_prev
                        inv_prev.loc[inv_exp_solar_prev.index, ("Solar", "capex")] += inv_exp_solar_prev


                L_selection = L_calc.loc[:, L_calc.columns.isin(investment_w_curr.index.astype(str))]
                L_selection.index = L_selection.index.astype(str)
                L_selection.columns = L_selection.columns.astype(str)

                L_prev = L_b.loc[:, L_b.columns.isin(investment_w_prev.index.astype(str))]
                L_prev.index = L_prev.index.astype(str)
                L_prev.columns = L_prev.columns.astype(str)

                # $x = (I-A)^{-1} F$, and so $\Delta x = L  \Delta F$

                capex_output_year = L_selection.dot(inv.loc[L_selection.columns, :]).drop(
                    "GFGD"
                )
                capex_output_year_prev = L_prev.dot(inv_prev.loc[L_prev.columns, :]).drop(
                    "GFGD"
                )

                capex_output_year -= capex_output_year_prev


                # Back to 2.1 Add batteries opex using capex data
                # annual change in generation spending
                # We want to add an additional column for batteries spending
                yr_b_usd = gen_scen_batt.set_index(["year", "scenario"]).loc[
                    (yr_b, global_scen)
                ]["annual.opex"]
                yr_e_usd = gen_scen_batt.set_index(["year", "scenario"]).loc[
                    (yr_e, global_scen)
                ]["annual.opex"]

                # average annual change in battery spending
                #annual_battery_delta = (yr_e_usd - yr_b_usd) / (yr_e - yr_b)

                # annual battery impact = Leontief (X) (battery investment weights * average annual change in spending)
                batt_out_year = (
                    L_calc.loc[:, L_calc.columns.isin(investment_w_curr.loc[:, ("Battery storage", "capex")].index.astype(str))]
                    .dot(investment_w_curr.loc[:, ("Battery storage", "capex")] * yr_e_usd)
                    .drop("GFGD")
                )
                batt_out_year_prev = (
                    L_b.loc[:, L_b.columns.isin(investment_w_prev.loc[:, ("Battery storage", "capex")].index.astype(str))]
                    .dot(investment_w_prev.loc[:, ("Battery storage", "capex")] * yr_b_usd)
                    .drop("GFGD")
                )

                # we measure in thousands, and add to opex_output_year
                opex_output_year["Battery storage"] = (batt_out_year - batt_out_year_prev) * 1000

                # 3 INDUSTRY OUTPUT TO OCCUPATION DEMAND
                # disaggregate all energy tech

                elec_names = [
                    "Hydroelectric power generation",
                    "Nuclear electric power generation",
                    "Solar electric power generation",
                    "Wind electric power generation",
                    "Geothermal electric power generation",
                    "Coal electric power generation",
                    "Biomass electric power generation",
                    "Gas electric power generation",
                ]
                elec_names_short = {
                    "Hydroelectric power generation": "Hydro",
                    "Nuclear electric power generation": "Nuclear",
                    "Solar electric power generation": "Solar",
                    "Wind electric power generation": "Wind",
                    "Geothermal electric power generation": "Geothermal",
                    "Coal electric power generation": "Coal",
                    "Gas electric power generation": "Natural gas",
                    "Biomass electric power generation": "Biomass",
                }

                # capex
                capex_output_year.drop("imports", inplace=True)
                capex_output_year.columns = capex_output_year.columns.droplevel(1)

                # opex
                opex_output_year.drop(["GFGD", "imports"], inplace=True)
                opex_output_year = opex_output_year.rename(elec_names_short, axis="columns")

            
                typ, spnding = ("capex", capex_output_year)
                for tech in spnding.columns.drop_duplicates():
                    spending_vector = spnding.loc[:, tech].copy()
                    if (tech == 'T&D') & (settings.SENS_TANDD_CONSTR_BREAKDOWN == True):
                        emp_tandd_constr = help_sens.tandd_constr_sens(spending_vector.loc['23'], emp_per_million_ouput_sum_2018, settings)
                        spending_vector.loc['23'] = 0

                    emp_invest[tech + "-" + typ] = (
                        emp_per_million_ouput_sum_2018.loc[:, spending_vector.index]
                        @ spending_vector
                    )

                    if (tech == 'T&D') & (settings.SENS_TANDD_CONSTR_BREAKDOWN == True):
                        emp_invest[tech + "-" + typ] += emp_tandd_constr

                    emp_invest[tech + "-" + typ] /= 1000000

                typ, spnding = ("opex", opex_output_year)
                for tech in spnding.columns.drop_duplicates():
                    spending_vector = spnding.loc[:, tech].copy()
                    emp_invest[tech + "-" + typ] = (
                        emp_per_million_ouput_sum_2018.loc[:, spending_vector.index]
                        @ spending_vector
                    )
                    emp_invest[tech + "-" + typ] /= 1000000

                # reset or initialize emp_tot
                emp_invest["emp_tot"] = 0
                emp_invest["emp_tot"] = emp_invest.drop(
                    [
                        "TOT_EMP",
                        "A_MEAN",
                        "TOT_EMP_minor",
                        "A_MEAN_minor",
                        "TOT_EMP_major",
                        "A_MEAN_major",
                        'OCC_TITLE', 'class_minor', 'class_major', 'OCC_CODE_minor',
                        'OCC_TITLE_minor', 'OCC_CODE_major', 'OCC_TITLE_major',
                        'classification', 'color'
                    ],
                    axis=1,
                ).sum(axis=1)

                # 4 SAVE DATA
                # save sensitivity analysis summary data
                new_year = str(yr_b) + "-" + str(yr_e)
                #new_jobs[new_year] = np.nan 
            
                # calculate each year's cumulative employment total
                # make the latest year's jobs a cummulative total
                if yr_b > min(yrs):
                    old_year = str(yr_b - 1) + "-" + str(yr_b)
                    new_jobs[new_year] = emp_invest['emp_tot'] + new_jobs[old_year]
                else: #else it just equals the first year emp_tot
                    new_jobs[new_year] = emp_invest['emp_tot']
           
                # occupations
                if settings.sub_run_id ==1: #only save the first sub-run if there are multiple
                    emp_invest.to_csv(
                        settings.FL_ann_data_out
                        + str(MODELRUN_ID)
                        + "_empl_"
                        + global_scen
                        + "_"
                        + str(yr_b)
                        + "-"
                        + str(yr_e)
                        + ".csv"
                    )
                if settings.sub_run_id > 1: #save other subruns with sun_run_id 
                    emp_invest.to_csv(
                        settings.FL_ann_data_out
                        + str(MODELRUN_ID)
                        + '_' + str(settings.sub_run_id)
                        + "_empl_"
                        + global_scen
                        + "_"
                        + str(yr_b)
                        + "-"
                        + str(yr_e)
                        + ".csv"
                    )


                ## Industries
                capex_ind_workers = (
                    capex_output_year.T * emp_per_million_ouput_sum_2018.sum() / 1000000
                )
                opex_ind_workers = (
                    opex_output_year.T * emp_per_million_ouput_sum_2018.sum() / 1000000
                )

                capex_ind_workers.index += "_capex"
                opex_ind_workers.index += "_opex"

                if settings.sub_run_id ==1: #only save the first sub-run if there are multiple
                    pd.concat([capex_ind_workers, opex_ind_workers]).T.to_csv(
                        settings.FL_ann_data_out
                        + str(MODELRUN_ID)
                        + "_ind_empl_"
                        + global_scen
                        + "_"
                        + str(yr_b)
                        + "-"
                        + str(yr_e)
                        + ".csv"
                    )
                    
                if settings.sub_run_id > 1: #save other subruns with sun_run_id 
                    pd.concat([capex_ind_workers, opex_ind_workers]).T.to_csv(
                        settings.FL_ann_data_out
                        + str(MODELRUN_ID)
                        + '_' + str(settings.sub_run_id)
                        + "_ind_empl_"
                        + global_scen
                        + "_"
                        + str(yr_b)
                        + "-"
                        + str(yr_e)
                        + ".csv"
                    )

            #calculate the maximum jobs added for each occupation and scenario from new_jobs output
            years = new_jobs.shape[1]
            #-	2034 net new jobs since 2020
            new_jobs['max_2034_jobs'] = new_jobs.iloc[:, 2:15].max(1)
            #-	2048 net new jobs since 2020
            new_jobs['max_2048_jobs'] = new_jobs.iloc[:, 2:years].max(1)
            #-49-9080	Maximum new wind turbine technicians any year since 2020
            #-47-2230	Maximum new solar PV installers any year since 2020
            #-47-5040	Maximum mining machine operators lost any year since 2020
            #-51-8010	Maximum power plant operators jobs lost any year since 2020
            key_jobs = new_jobs.loc[["49-9080","47-2230","47-5040","51-8010", "49-9051"]]

            # save output summary file for every subrun 
            new_jobs.to_csv(settings.FL_data_out + str(MODELRUN_ID) + '_' + str(settings.sub_run_id) + "_new_jobs_"  
                 + "_" + global_scen + ".csv")
            key_jobs.to_csv(settings.FL_data_out + str(MODELRUN_ID) + '_' + str(settings.sub_run_id) + "_key_jobs_"  
                 + "_" + global_scen + ".csv")

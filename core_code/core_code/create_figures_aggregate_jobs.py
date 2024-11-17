import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import settings
import help_sens
from itertools import cycle, islice

#### JORIS - these are all in settings
# 0 SET PARAMETERS
#analysis = "../../../../analysis/"
#data = "../../../../data/"
#data_out = "../../results_sens/annual_change_sens"
#fig_folder = "../../figs/"
#version = "3_oct_2022"

def create_figs(MODELRUN_ID, settings):

#get the years smoothing according to the sensitivity analysis that is being run
    # if not smoothing than smoothing years = 1
    if settings.SENS_ANALYSIS_SMOOTHING == True:
        years_smoothing = settings.YEARS_SMOOTHING_SENS
        smoothing = True
    else:
        years_smoothing = 1
        smoothing = False

    ylim = (-250000, 600000)  # for the cumulative area plots
    if smoothing == False:
        yrs = range(2021, 2049)
        local_data_out = settings.FL_ann_data_out + "_no_smoothing/"
    else:
        local_data_out = settings.FL_ann_data_out

    scens = ["95% by 2050", "95% by 2035", "Reference"]
    fast_scen = "95% by 2035"
    ref_scen = "Reference"


    # FIGURE 1
    for scen in scens:
        df_full = pd.DataFrame()
        for yr in range(2021, 2048, 1):
            df = pd.read_csv(
                local_data_out
                + str(MODELRUN_ID)
                + "_ind_empl_"
                + scen
                + "_"
                + str(yr)
                + "-"
                + str(yr + 1)
                + ".csv",
                index_col=0,
            )
            df["year"] = yr
            df_full = pd.concat([df_full, df])

        df_full = df_full.reset_index().set_index(["index", "year"])
        df_full.groupby(df_full.index.get_level_values("year")).sum().cumsum().T.plot.bar(
            figsize=(15, 5)
        ).legend(loc="upper right", bbox_to_anchor=(1.07, 1.0), prop={"size": 7})
        plt.tight_layout();

        # Only create figure if settings.CREATE_FIGURES == True:
        if settings.CREATE_FIGURES == True:
            plt.savefig(
                settings.FL_fig_folder
                + str(MODELRUN_ID)
                + "_jobs_per_tech_opex_capex_scen_"
                + scen
                + "_smooth_"
                + str(smoothing)
                + "_2020-2050.png"
            )


    # FIGURE 2
    # for Matt, until 2035
    df_full = pd.DataFrame()
    for yr in range(2022, 2035, 1):
        df = pd.read_csv(
            local_data_out
            + str(MODELRUN_ID)
            + "_ind_empl_"
            + fast_scen
            + "_"
            + str(yr)
            + "-"
            + str(yr + 1)
            + ".csv",
            index_col=0,
        )
        df["year"] = yr
        df_full = pd.concat([df_full, df])

    df_full = df_full.reset_index().set_index(["index", "year"])
    df_full.groupby(df_full.index.get_level_values("year")).sum().cumsum().T.plot.bar(
        figsize=(15, 5)
    ).legend(loc="upper right", bbox_to_anchor=(1.07, 1.0), prop={"size": 7})

    plt.title(scen)
    plt.tight_layout();

    # Only create figure if settings.CREATE_FIGURES == True:
    if settings.CREATE_FIGURES == True:
        plt.savefig(
            settings.FL_fig_folder
            + str(MODELRUN_ID)
            + "_jobs_per_tech_opex_capex_scen_"
            + scen
            + "_smooth_"
            + str(smoothing)
            + "_2020-2035.png"
        )


    # FIGURE 3
    for scen in scens:
        df_full = pd.DataFrame()
        for yr in range(2021, 2048, 1):
            df = pd.read_csv(
                local_data_out
                + str(MODELRUN_ID)
                + "_ind_empl_"
                + scen
                + "_"
                + str(yr)
                + "-"
                + str(yr + 1)
                + ".csv",
                index_col=0,
            )
            df["year"] = yr
            df_full = pd.concat([df_full, df])

        df_full = df_full.reset_index().set_index(["index", "year"])
        df_try = df_full.copy()
        df_try.index = df_try.index.droplevel(0)
        df_try = df_try.groupby(df_try.index).sum()
        df_try.loc[2020, :] = 0
        df_try.sort_index(inplace=True)

        my_colors = list(
            islice(
                cycle(
                    [
                        "g",
                        "yellow",
                        "grey",
                        "k",
                        "lime",
                        "peru",
                        "b",
                        "indigo",
                        "darkorange",
                        "palegreen",
                        "k",
                        "lightyellow",
                        "slategrey",
                        "peachpuff",
                        "dodgerblue",
                        "fuchsia",
                        "mediumseagreen",
                        "gainsboro",
                        "orange",
                    ]
                ),
                None,
                2 * len(df_try.columns),
            )
        )

        pd.concat(
            [df_try.cumsum().clip(lower=0), df_try.cumsum().clip(upper=0)], axis=1
        ).plot.area(color=my_colors, figsize=(6, 6), ylim=ylim, legend="reverse").legend(
            loc="center left", bbox_to_anchor=(1.0, 0.5)
        )

        plt.ylabel("Cumulative new jobs in electricity supply chain")
        plt.title(scen)
        plt.tight_layout();
        
        # Only create figure if settings.CREATE_FIGURES == True:
        if settings.CREATE_FIGURES == True:
            plt.savefig(
                settings.FL_fig_folder
                + str(MODELRUN_ID)
                + "_jobs_cumul_stack_per_tech_opex_capex_scen_"
                + scen
                + "_smooth_"
                + str(smoothing)
                + "_2020-2050.png"
            )

    ind_agg = pd.read_csv(
        settings.FL_data + "Data_out_names/ind_sum_aggr.csv",
        index_col=0,
    )["label aggr"]


    # FIGURE 4
    for scen in scens:
        df_full = pd.DataFrame()
        for yr in range(2021, 2048, 1):
            df = pd.read_csv(
                local_data_out
                + str(MODELRUN_ID)
                + "_ind_empl_"
                + scen
                + "_"
                + str(yr)
                + "-"
                + str(yr + 1)
                + ".csv",
                index_col=0,
            )
            df["year"] = yr
            df_full = pd.concat([df_full, df])

        df_full = df_full.reset_index().set_index(["index", "year"])
        df_try = df_full.unstack().T.copy()
        df_try.index = df_try.index.droplevel(0)
        df_try = df_try.groupby(df_try.index).sum()
        df_try.loc[2020, :] = 0
        df_try.sort_index(inplace=True)

        utils = set(df_try) - set(ind_agg.index)

        utils = dict(zip(utils, ["Utilities" for x in utils]))

        ind_agg = pd.concat([ind_agg, pd.Series(utils, dtype=str)])

        indcat = np.unique(ind_agg)
        my_colors = list(
            islice(
                cycle(
                    [
                        "lime",
                        "violet",
                        "y",
                        "darkred",
                        "aquamarine",
                        "b",
                        "r",
                        "pink",
                        "peru",
                        "k",
                        "grey",
                        "slategrey",
                        "midnightblue",
                        "salmon",
                        "goldenrod",
                        "teal",
                        "purple",
                        "deeppink",
                        "g",
                        "m",
                        "c",
                    ]
                ),
                None,
                len(indcat),
            )
        )
        colordict = dict(zip(indcat, my_colors))
        colorlist = [colordict[ind_agg.to_dict()[ind]] for ind in df_try.columns]

        #     my_colors = list(islice(cycle(['b', 'r', 'g', 'y', 'k', 'm', 'c', 'lime',
        #                                    'slategrey', 'midnightblue', 'salmon', 'goldenrod',
        #                                    'teal', 'purple', 'violet', 'pink', 'peru', 'grey']), None, 2*len(df_try.columns)))

        pd.concat(
            [df_try.cumsum().clip(lower=0), df_try.cumsum().clip(upper=0)], axis=1
        ).plot.area(
            color=colorlist, figsize=(14,7), ylim=ylim, legend="reverse"
        )  # .\
        # legend(loc='center left', bbox_to_anchor=(1.0, 0.5))

        markers = [
            plt.Line2D([0, 0], [0, 0], color=color, marker="o", linestyle="")
            for color in colordict.values()
        ]
        plt.legend(
            markers,
            colordict.keys(),
            numpoints=1,
            loc="center left",
            bbox_to_anchor=(1.0, 0.5),
        )

        plt.ylabel("Cumulative worker demand in the electricity supply chain")
        plt.title(scen)
        plt.tight_layout();

        # Only create figure if settings.CREATE_FIGURES == True:
        if settings.CREATE_FIGURES == True:
            plt.savefig(
                settings.FL_fig_folder
                + str(MODELRUN_ID)
                + "_jobs_cumul_stack_per_ind_opex_capex_scen_"
                + scen
                + "_smooth_"
                + str(smoothing)
                + "_2020-2050.png"
            )


    # FIGURE 5
    major_occ = pd.read_csv(
        local_data_out + str(MODELRUN_ID) +  "_empl_95% by 2035_2021-2022.csv",
        index_col=0,
    )["OCC_TITLE_major"]

    cols = [
        "Wind-capex",
        "Solar-capex",
        "Natural gas-capex",
        "Coal-capex",
        "Biomass-capex",
        "Geothermal-capex",
        "Hydro-capex",
        "Battery storage-capex",
        "T&D-capex",
        "Biomass-opex",
        "Coal-opex",
        "Solar-opex",
        "Natural gas-opex",
        "Geothermal-opex",
        "Hydro-opex",
        "Nuclear-opex",
        "Wind-opex",
        "Transmission and distribution-opex",
        "Battery storage-opex",
    ]

    for scen in scens:
        df_full = pd.DataFrame()
        for yr in range(2021, 2048, 1):
            df = pd.read_csv(
                local_data_out
                + str(MODELRUN_ID)
                + "_empl_"
                + scen
                + "_"
                + str(yr)
                + "-"
                + str(yr + 1)
                + ".csv",
                index_col=0,
            )
            df["year"] = yr
            df_full = pd.concat([df_full, df])

        df_full = df_full.reset_index().set_index(["index", "year"])
        df_try = df_full.loc[:, cols].unstack().T.copy()
        df_try.index = df_try.index.droplevel(0)
        df_try = df_try.groupby(df_try.index).sum()
        df_try.loc[2020, :] = 0
        df_try.sort_index(inplace=True)

        occcat = np.unique(major_occ.str[:-12])
        my_colors = list(
            islice(
                cycle(
                    [
                        "b",
                        "r",
                        "g",
                        "y",
                        "k",
                        "m",
                        "c",
                        "lime",
                        "slategrey",
                        "midnightblue",
                        "salmon",
                        "goldenrod",
                        "teal",
                        "purple",
                        "violet",
                        "pink",
                        "peru",
                        "grey",
                        "deeppink",
                        "olive",
                        "darkred",
                        "aquamarine",
                    ]
                ),
                None,
                len(occcat),
            )
        )
        colordict = dict(zip(occcat, my_colors))
        colorlist = [colordict[occ] for occ in major_occ.str[:-12]]

        pd.concat(
            [df_try.cumsum().clip(lower=0), df_try.cumsum().clip(upper=0)], axis=1
        ).plot.area(
            color=colorlist, figsize=(11,7), ylim=ylim, legend="reverse"
        )  # .\
        # legend(loc='center left', bbox_to_anchor=(1.0, 0.5))

        markers = [
            plt.Line2D([0, 0], [0, 0], color=color, marker="o", linestyle="")
            for color in colordict.values()
        ]
        plt.legend(
            markers,
            colordict.keys(),
            numpoints=1,
            loc="center left",
            bbox_to_anchor=(1.0, 0.5),
        )

        plt.ylabel("Cumulative new jobs in electricity supply chain")
        plt.title(scen)
        plt.tight_layout();

        # Only create figure if settings.CREATE_FIGURES == True:
        if settings.CREATE_FIGURES == True:
            plt.savefig(
                settings.FL_fig_folder
                + str(MODELRUN_ID)
                + "_jobs_cumul_stack_per_occ_opex_capex_scen_"
                + scen
                + "_smooth_"
                + str(smoothing)
                + "_2020-2050.png"
            )
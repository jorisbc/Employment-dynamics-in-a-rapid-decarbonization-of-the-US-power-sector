energy_scenario_code
- Standalone R project to calculate the cost from NREL's energy scenarios

data
- all input data for energy-scenario-cost-to-occupational-frictions analysis, either from source (BLS, BEA) or from earlier analysis

results:
- data_out: output data from analysis, both intermediary and final
- figs: main figures


Order of execution:
1) preprocessing notebooks
- Jupyter notebooks to transform input data for main routines

2) core_code
- Routines to calculate occupation demand per year

3) postprocessing notebooks
- Jupyter notebooks to make the figures in the paper

4) suplemental materials notebooks
- Figures and analysis for supplemental materials

See further READMEs in the different folders
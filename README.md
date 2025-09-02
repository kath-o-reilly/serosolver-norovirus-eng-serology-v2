# serosolver-norovirus-english-serology

Analysis and data for serosolver analysis of english data

The analysis is in a series of sections.

NOTE: run files in R while in the `serosolver-norovirus-english-serology` Project

# 1. Check the serosolver vingette works

The norovirus application is very similar to the second example (cross-sectional 'flu data'), so check this example works ok.

In the example there are serology data from 1359 individuals ('titre_dat') that are analysed. This data are more comprehensive than norovirus: larger sample size and more 'flu variants.

Corresponding file: `serosolver_example2_vignette.R`

Note that `serosolver_example2_norovirus.R` is present, but refers to data found elsewhere on computer.

# 2. Data curation

File `plot_figure_1.R` should load in the norovirus variant data and serology samples that make Figure 1. (*Need to check*).

File `make_plots_from_serology_data.R` should load data and make dataframes for loading into serosolver (*Need to check*).

# 3. Creating AC maps

File `make_ac_maps.R` should load in data and make the AC maps with right colours (*Need to check*)

# 4. Run serosolver on the norovirus data

Sevral analyses are carried out;
- Blockade data and Debbink cartography

Directly generates parameter estimates, figure 3c (attack rates)

- Blockade data and Kendra cartography
- Blockade and avidity data applied to Debbink cartography `2008-2012_verified_IC50_23AUG2020.csv`

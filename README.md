# traits-and-ecosystem-properties
R scripts associated to PLANT TRAITS ALONE ARE POOR PREDICTORS OF ECOSYSTEM PROPERTIES AND LONG-TERM ECOSYSTEM FUNCTIONING manuscript

main analysis.r: data files (see Pangea) are read. Based on community and trait data, FD and FI indicators are calculated for each plot in each year. Then, all eacosystem properties are analysed. Proportion of variance explained by FI and FD is calculated for each function. Then, overlap in significant predictors for each combination of ecosystem functions is calculated. The same is done for PCA values based on FD/FI. Latly, only the 6 traits most frequently assessed in other studies are used to explain ecosystem properties.

one random trait analysed.r: 41 analyses are done, in which in each of them another of the 41 traits is analysed to explain ecosystem properties. R2 values of models are extracted

random sets of traits analysed.r: A random subset of the 41 traits (between 2 and 39 traits) are analysed as predictors of ecosystem functioning. 100 random subsets are taken. For each number of traits, R2 values of models are extracted.

forty random trait analysed.r: 41 analyses are done, in which in each of them another subset of 40 of the 41 traits is analysed to explain ecosystem properties. R2 values of models are extracted.

final plot.r: R2 values extracted from the above three R files are extracted, and correlated to the number of traits that is analysed.

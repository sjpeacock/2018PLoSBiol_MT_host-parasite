### Supporting code for 
## [Empirical evidence that metabolic theory describes the temperature dependency of within-host parasite dynamics](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2004608 "Link to article")  
### Devin Kirk*, Natalie Jones, Stephanie Peacock**, Jessica Phillips, Peter K. Molnar, Martin Krkosek, and Pepijn Luijckx
 '*' corresponding author: devin.kirk@mail.utoronto.ca

** questions regarding this code should be directed to stephanie.j.peacock@gmail.com

**Last updated**: Feb 11, 2018

**Status**: Code should be uploaded, but workspaces with final output are not provided because they are too big. Email me if you want to be able to reproduce the figures etc. and I can share them with you some other way.

#### This repository contains:

1. **DT**: folder that contains the code to fit the Discrete Temperature (DT) model where each temperature treatment is fit to the data independently. 
 - `DT_model_fitting.R`: Loads the `mortality_data.csv` and fits the DT model to each temperature in a loop, also looping over 15 'clones' of the data. Runs 10 independent chains for each fit, on parallel processors (so requires 10 processors to run as written). See Supporting Information for details. Requires *JAGS* to be installed and the package `dclone`.
 - `DT_analysing_output.R`: Loads a workspace file saved from `DT_model_fitting.R` to look at the chains for convergence, assess estimability, and output final parameter values, etc.
 - `DT_priors.csv`: File that contains the mean and sd of the log-normal prior distributions for DT model parameters. See also Table S2 of Supporting Information.
 - `DT_ParameterEstimates_2017Apr4_lognormal.xlsx`: Excel file containing the output of the final model fits.

2. **MTE**: folder that contains the code to fit the Metabolic Theory of Ecology (MTE) model where the model parameters at each temperature are given by a metabolic function, whose hyperparameters are estimated here.
- `MTE_model_fitting.R`: Loads the `mortality_data.csv` and fits the MTE model to each temperature in a loop, also looping over 25 'clones' of the data. Runs 10 independent chains for each fit, on parallel processors (so requires 10 processors to run as written). See Supporting Information for details. Requires *JAGS* to be installed and the package `dclone`.
- `MTE_analysing_output.R`: Loads a workspace file saved from `DT_model_fitting.R` to look at the chains for convergence, assess estimability, and output final parameter values, etc.
- `MTE_priors.csv`: File that contains the mean and sd of the log-normal prior distributions for DT model parameters. See also Table S3 of Supporting Information.
 - `CI_calc.R`: Monte Carlo simulation to produce confidence intervals on MTE predictions for model parameters, shown in Fig 3 of the main text. This requires the R workspace with output from `MTE_model_fitting.R`, which I can provide on request (too large to post to GitHub).

3. **`mortality_data.csv`**: Data from the experiments used for model fitting. Each row describes an individual host Daphnia. Variables include:
 - `temp`: the approximate temperature in degrees Celcius of the treatment.
- `treatment`: dummy variable for whether the individual host was exposed (1) or unexposed (0) to parasites.
- `infection_status`: whether the host was exposed (Infected) or unexposed (Uninfected) to parasites. Directly corresponds to treatment variable.
- `ID  `: unqiue identifier for individual host.
- `total_spvs`: number of spore clusters in a host at death. Eight individuals were known to be infected (see variable secondary below) but could not be dissected, and thus have a total_spvs of `NA`.
- `death`: day of the experiment when the host was found dead. We assumed that the host died at some point between this day and the previous day, since hosts were examined daily for signs of vitality.
- `secondary`: Whether a secondary infection of uninfected juveniles was successful, indicating that the deceased individual was infected even though we were unable to successfully dissect it and the parasite load at death was unknown (total_spvs = `NA`). See Methods. 

4. **`models.R`**: JAGS code for the DT and MTE models.

5. **`figures.R`**: R code to reproduce the figures in the paper. To run this code, you need to load R workspaces with the model results, which were too big to post to GitHub (>100 MB). The code for the figures is provided for transparency, but if you actually want to run it then email me and I can send the workspaces.

6. **`plotting_functions.R`**: some functions used in the analysing_output files and in the figures.

## Note on notation:
There is a discrepency between the manuscript and the code in notation, due to a last-minute change in notation in the published manuscript. The parameter theta describing the carrying capacity of the parasite popultion in the manuscript is written as K throughout the code. Also, the DT model is written as TD in much of the code. 

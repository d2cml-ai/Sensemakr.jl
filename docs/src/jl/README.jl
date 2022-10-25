# # Sensemakr
#
# [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://d2cml-ai.github.io/Sensemakr.jl/stable/)
# [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://d2cml-ai.github.io/Sensemakr.jl/dev/)
# [![Build Status](https://github.com/d2cml-ai/Sensemakr.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/d2cml-ai/Sensemakr.jl/actions/workflows/CI.yml?query=branch%3Amaster)
# [![Coverage](https://codecov.io/gh/d2cml-ai/Sensemakr.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/d2cml-ai/Sensemakr.jl)
#
# sensemakr for Julia (`Sensemakr`) implements a suite of sensitivity analysis tools that extends the traditional omitted variable bias framework and makes it easier to understand the impact of omitted variables in regression models, as discussed in Cinelli, C. and Hazlett, C. (2020) “Making Sense of Sensitivity: Extending Omitted Variable Bias.” Journal of the Royal Statistical Society, Series B (Statistical Methodology).
#
# ## Installation
#
# You will need Julia version 1.7.0 or higher to install this package
#
# The latest version can be downloaded by running:
#
# ````julia
# import Pkg; Pkg.add(url = "https://github.com/d2cml-ai/Sensemakr.jl");
# ````
#
# The latest stable version can be downloaded by running:
#
# ````julia
# import Pkg; Pkg.add("Sensemakr");
# ````
#
# ## Example Usage

## loads Sensemakr
using Sensemakr

## loads DataFrames and GLM
using DataFrames, GLM

## creates dataframe with darfur data
darfur = load_darfur();

## runs OLS
form = @formula(peacefactor ~ directlyharmed + age + farmer_dar + herder_dar + pastvoted + hhsize_darfur + female + village);
fitted_model = lm(form, darfur);

## creates sensemakr object
darfur_sense = sensemakr(fitted_model, "directlyharmed", benchmark_covariates = "female", kd = [1, 2, 3]);

## summary of sensitivity analysis
summary(darfur_sense)

#

## contour plot for estimate
plot(darfur_sense)

# ![Figure_1](docs/src/images/Figure_1.png)

## contour plot for t-value
plot(darfur_sense, sensitivity_of = "t-value")

# ![Figure_2](docs/src/images/Figure_2.png)

## extreme scenario plot
plot(darfur_sense, plot_type = "extreme")

# ![Figure_3](docs/src/images/Figure_3.png)


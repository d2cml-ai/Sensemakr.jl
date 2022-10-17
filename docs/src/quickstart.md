# Quickstart

## Introduction

In this quickstart, we demonstrate the basic usage of the package by reproducing Section 5 of Cinelli and Hazlett (2020), which estimates the effects of exposure to violence on attitudes towards peace, in Darfur. Throughout this manual, we mainly focus on the code. For detailed explanations, please refer to the original paper.

## Violence in Darfur

In 2003 and 2004, the Darfurian government orchestrated a horrific campaign of violence against civilians, killing an estimated two hundred thousand people. In this application, we are interested in answering the following question: did being directly exposed to harm make individuals more “angry,” and thus more likely to ask for revenge, or did it make them more "weary," and thus more likely to ask for peace?

`Sensemakr` comes with the Darfur dataset, which can be loaded with the command `data.load_darfur()`. More details about the data can be found in the documentation, as well as in Hazlett (2019) and Cinelli and Hazlett (2020).

```julia
# imports Sensemakr
julia> using Sensemakr

# loads darfur data
julia> using DataFrames

julia> darfur = load_darfur();

julia> first(darfur, 5)
5×14 DataFrame
 Row │ wouldvote  peacefactor  peace_formerenemies  peace_jjindiv  peace_jjtribes  gos_soldier_execute  directlyharmed ⋯
     │ Int64      Float64      Int64                Int64          Int64           Int64                Int64          ⋯
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │         0     1.0                         1              1               1                    0               0 ⋯
   2 │         0     0.706831                    0              1               1                    0               0
   3 │         1     0.0                         0              0               0                    1               0
   4 │         1     0.495178                    1              0               0                    0               1
   5 │         0     0.0                         0              0               0                    1               1 ⋯
```

A natural approach for such problem is to run the following linear regression model, where we regress `peacefactor` on `directlyharmed`, further adjusting for `village`, `female` as well as other covariates. Here we run this regression using `GLM`.

```julia
julia> using GLM

julia> form = @formula(peacefactor ~ directlyharmed + age + farmer_dar + herder_dar + pastvoted + hhsize_darfur + female + village);

julia> fitted_model = lm(form, darfur);
```

The above regression results in the following estimate and standard errors for the coefficient of `directlyharmed`:

```julia
julia> first(coeftable(fitted_model), 2)[2]
(Name = "directlyharmed", var"Coef." = 0.09731581928495726, var"Std. Error" = 0.023256537809811153, t = 4.184449984808271, var"Pr(>|t|)" = 3.182339959789431e-5, var"Lower 95%" = 0.05166327477010026, var"Upper 95%" = 0.14296836379981426)
```

According to this model, those individual who were directly exposed to harm became on average more “pro-peace,” not less.

## Sensitivity Analysis

The causal interpretation of the previous estimate, however, relies on the assumption that `village` and `gender` are sufficient for control of confounding—in other words, it requires the assumption of no unobserved confounders. What if is wrong? How strong would these unobserved variables have to be in order to change the original research conclusions? The goal of `Sensemakr` is precisely that, i.e, to make it easier to understand the impact that omitted variables would have on a regression result.

The main function of the package is the function `sensemakr`. The `sensemakr` function defines an object of type `sensemakr`, performing the most commonly required sensitivity analyses which can then be further explored with the `summary` and `plot` methods of the object. This function is mainly a convenience wrapper for other sensitivity functions defined in the package, which can also be called directly, as we detail later in the documentation.

In the code chunk below, we apply the function `sensemakr` to our OLS model from statsmodel.

```julia
# creates a sensemakr object
julia> darfur_sense = sensemakr(fitted_model, "directlyharmed", benchmark_covariates = "female", kd = [1, 2, 3], ky = [1, 2, 3], q = 1.0, alpha = 0.05, reduce = true);
```

The main arguments of the call are:

**model**: the OLS model with the outcome regression. In our case, `darfur_model`.

**treatment**: the name of the treatment variable. In our case, "directlyharmed".

**benchmark_covariates**: the names of covariates that will be used to bound the plausible strength of the unobserved confounders. Here, we put "female", which is arguably one of the main determinants of exposure to violence, and also a strong determinant of attitudes towards peace.

**kd** and **ky**: these arguments parameterize how many times stronger the confounder is related to the treatment (kd) and to the outcome (ky) in comparison to the observed benchmark covariates (in this case, female). In our example, setting `kd = [1, 2, 3]` and ky = `[1, 2, 3]` means we want to investigate the maximum strength of a confounder once, twice, or three times as strong as female (in explaining treatment and outcome variation). If only `kd` is given, `ky` will be set equal to `kd`.

**q**: fraction of the effect estimate that would have to be explained away to be problematic. Setting q = 1, means that a reduction of 100% of the current effect estimate, that is, a true effect of zero, would be deemed problematic. The default is 1.

**alpha**: significance level of interest for statistical inference. The default is 0.05.

**reduce**: should we consider confounders acting towards increasing or reducing the absolute value of the estimate? The default is `reduce = true`, which means we are considering confounders that pull the estimate towards (or through) zero.

Using the default arguments, one can simplify the previous call to:

```julia
julia> darfur_sense = sensemakr(fitted_model, "directlyharmed", benchmark_covariates = "female", kd = [1, 2, 3]);
```

Once we run `sensemakr`, we can now explore the sensitivity analysis results.


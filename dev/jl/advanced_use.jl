# # Advanced usage
#
# The standard functionality of the `Sensemakr`, demonstrated in the previous section, will suﬃce for most users, most of the time. When needed, more ﬂexibility can be obtained by employing the sensitivity functions of the package directly, such as: (i) functions for sensitivity plots (`ovb_contour_plot`, `ovb_extreme_plot`); (ii) functions for computing bias-adjusted estimates and t-values (`adjusted_estimate`, `adjusted_t`); (iii) functions for computing the robustness value and partial R2 directly (`robustness_value`, `partial_r2`); and, (iv) functions for bounding the strength of unobserved confounders (`ovb_bounds`), among others. These functions are automatically imported when loading the package.
# 
# In this page we demonstrate how to use these functions in two examples.
#
# ## The risks of informal benchmarking
#
# Informal “benchmarking” procedures have been widely suggested in the sensitivity analysis literature as a means to aid interpretation. It intends to describe how an unobserved confounder $Z$ “not unlike” some observed covariate $X_j$ would alter the results of a study. [Cinelli and Hazlett (2020)](https://doi.org/10.1111/rssb.12348) show why these proposals may lead users to erroneous conclusions, and offer formal bounds on the bias that could be produced by unobserved confounding “as strong” as certain observed covariates.
#
# Here we use `Sensemakr` to replicate the example in Section 6.1 of [Cinelli and Hazlett (2020)](https://doi.org/10.1111/rssb.12348). This example provides a useful tutorial on how users can construct their own sensitivity contour plots with customized bounds, beyond what is offered by default on the package.
#
# ### Simulating the data
#
# We begin by simulating the data generating process used in our example. Consider a treatment variable $D$, an outcome variable $Y$, one observed confounder $X$, and one unobserved confounder $Z$. Again, all disturbance variables $U$ are standardized mutually independent gaussians, and note that, in reality, the treatment $D$ has no causal effect on the outcome $Y$.
#
# $$ Z = U_z \\
# X = U_x \\
# D = X + Z + U_d \\
# Y = X + Z + U_y $$
#
# Also note that, in this model: (i) the unobserved confounder `Z` is independent of `X`; and, (ii) the unobserved confounder `Z` is exactly like `X` in terms of its strength of association with the treatment and the outcome. The code below creates a sample of size 100 of this data generating process. We make sure to create residuals that are standardized and orthogonal so that all properties that we describe here will hold exactly even in this finite sample.

## loads GLM, DataFrames, and Statistics
using GLM, DataFrames, Statistics

## defines function to simulate orthogonal residuals
function resid_maker(n, df)
    N = randn(n)
    form = Term(:N) ~ sum(term.(propertynames(df)))
    df[:, "N"] = N
    model = lm(form, df)
    e = residuals(model)
    e = (e .- mean(e)) / std(e)
    return e
end;

## simulates data
n = 100;
X = randn(n);
X = (X .- mean(X)) / std(X);
Z = resid_maker(n, DataFrame(x = X));
D = X + Z + resid_maker(n, DataFrame(x = X, z = Z));
Y = X + Z + resid_maker(n, DataFrame(x = X, z = Z, d = D));

# ### Fitting the model
# In this example, the investigator does not observe the confounder $Z$. Therefore, she is forced to fit the restricted linear model $Y \sim D + X$, resulting in the following estimated values

## creates DataFrame with the data
df = DataFrame(Y = Y, X = X, Z = Z, D = D);

## runs OLS
model_ydx = lm(@formula(Y ~ D + X), df);

## print OLS table
println(coeftable(model_ydx))

## print R2
println(r2(model_ydx))

# Note we obtain a large and statistically significant coefficient estimate of the effect of $X$ on $Y$ ($0.5$). However, we know that the variable $Z$ is not observed, and there is the fear that this estimated effect is in fact due to the bias caused by $Z$. On the other hand, let us suppose the investigator correctly knows that: (i) Z and X have the same strength of association with $D$ and $Y$; and, (ii) $Z$ is independent of $X$. Can we leverage this information to understand how much bias a confounder $Z$ “not unlike” $X$ could cause?
#
# ### Informal benchmarks
#
# Computing the bias due to the omission of $Z$ requires two sensitivity parameters: its partial $R^2$ with the treatment $D$ and its partial $R^2$ with the outcome $Y$. How could we go about computing the bias that a confounder $Z$ “not unlike X” would cause?
#
# Intuitively, it seems that we could take as reference the observed partial $R^2$ of $X$ with $D$ and $Y$, and use those as the plausible values for the sensitivity parameters. So let us now compute those observed partial $R^2$ using the `partial_r2()`. For the partial $R^2$ of $X$ with the treatment, we also need to fit a treatment regression $D \sim X$ first.

## load Sensemakr
using Sensemakr

## fits treatment regression
model_dx = lm(@formula(D ~ X), df);

## computes partial R2 of X with D and Y
r2dx = partial_r2(model = model_dx, covariates = "X");
r2yx_d = partial_r2(model = model_ydx, covariates = "X");

# Once both partial $R^2$ are computed, we can determine the implied adjusted estimate due to an unobserved confounder $Z$ using the `adjusted_estimate()`.

## computes adjusted estimate
informal_adjusted_estimate = adjusted_estimate(r2dx, r2yx_d, model = model_ydx, treatment = "D");

# We can now plot the sensitivity contours with `ovb_contour_plot()` and add our informal benchmark with `add_bound_to_contour()`.

## draws sensitivity contours
ovb_contour_plot(model_ydx, "D", lim = 0.6)

## adds informal benchmark
add_bound_to_contour(r2dx, r2yx_d, bound_value = informal_adjusted_estimate, bound_label = ["Informal benchmark"], label_bump_x = 0.6/30, label_bump_y = 0.6/30)

# ![Figure_4](images/Figure_4.png)
#
# As we can see, the results of the informal benchmark are different from what we expected. The informal benchmark point is still far away from zero, and this would lead an investigator to incorrectly conclude that an unobserved confounder $Z$ “not unlike $X$” is not sufficient to explain away the observed effect. Moreover, this incorrect conclusion occurs despite correctly assuming that: (i) $Z$ and $X$ have the same strength of association with $D$ and $Y$; and, (ii) $Z$ is independent of $X$. See Section 6.1 of [Cinelli and Hazlett (2020)](https://doi.org/10.1111/rssb.12348) for details of why this happens.
#
# ### Formal Bounds
#
# We now show how to compute formal bounds using the function `ovb_bounds()`.

## computes formal bounds
formal_bound = ovb_bounds(model_ydx, "D", benchmark_covariates = "X", kd = 1, ky = 1);

# In this function we specify the linear model being used (`model_ydx`), the treatment of interest ($D$), the observed variable used for benchmarking ($X$), and how stronger $Z$ is in explaining treatment (`kd`) and outcome (`ky`) variation, as compared to the benchmark variable $X$. We can now plot the proper bound against the informal benchmark.

## contour plot contrasting informal vs formal bounds
ovb_contour_plot(model_ydx, "D", lim = 0.6)
add_bound_to_contour(r2dx, r2yx_d, bound_value = informal_adjusted_estimate, bound_label = ["Informal benchmark"], label_bump_x = 0.6/30, label_bump_y = 0.6/30)
add_bound_to_contour(formal_bound, bound_label = ["Proper bound"], label_bump_x = 0.6/30, label_bump_y = 0.6/30)

# ![Figure_5](images/Figure_5.png)
#
# Note that, using the formal bounds, the researcher now reaches the correct conclusion that, an unobserved confounder $Z$ similar to $X$ is strong enough to explain away all the observed association.
#
# ## Manually providing input data
#
# We now demonstrate how to replicate the sensitivity analysis of the Darfur example (see Quickstart) without access to the microdata, using only usual statistics found in the regression tables. All results of the Darfur example, except the formal benchmark bounds as we explaim below, can be recoverded by simply providing: (i) the point estimate of `directlyharmed` ($0.097$); (ii) its estimated standard error ($0.023$); and, (iii) the degrees of freedom of the regression ($783$). This type of “manual” analysis can be useful when reading or reviewing papers that did not perform sensitivity.

## computes the robustness value
robustness_value(t_statistic = 0.097 / 0.023, dof = 783)

#

## computes the partial R2
partial_r2(t_statistic = 0.097/0.023, dof = 783)

#

## computes what the adjusted estimate would have been
## after adjusting for confounders with strength r2dz_x = .1, r2yz_dx = 0.15
adjusted_estimate(0.1, 0.15, estimate = 0.097, se = 0.023, dof = 783)

#

## computes what the adjusted t-value would have been
## after adjusting for confounders with strength r2dz_x = .1, r2yz_dx = 0.15
adjusted_t(0.1, 0.15, estimate = 0.097, se = 0.023, dof = 783) 

#

## contour plot of the point estimate
ovb_contour_plot(0.097, 0.023, 783)

# ![Figure_6](images/Figure_6.png)

## contour plot of the t-value
ovb_contour_plot(0.097, 0.023, 783, sensitivity_of = "t-value")

# ![Figure_7](images/Figure_7.png)

## extreme plot
ovb_extreme_plot(0.097, 0.023, 783)

# ![Figure_8](images/Figure_8.png)
#
# Finally, we show how users can compute bounds on the strength of confounding using only summary statistics, if the paper also provides a treatment regression table, i.e., a regression of the treatment on the observed covariates. Such regressions are sometimes shown in published works as part of efforts to describe the “determinants” of the treatment, or as “balance tests” in which the investigator assesses whether observed covariates predict treatment assigment. For the Darfur example, the regression coefficient and standard errors (in parenthesis) of `female` for the outcome and treatment equation are -0.097 (0.036) and -0.232 (0.024). With these, we can now compute the formal bounds as shown below.

## computes partial R2 based on t-value of female
r2yxj_dx = partial_r2(t_statistic = -0.232 / 0.024, dof = 783);
r2dxj_x  = partial_r2(t_statistic = -0.097 / 0.036, dof = 783);

## computes bounds on the strength of Z if it were 1, 2 or 3 times as strong as female
bounds   = ovb_partial_r2_bound(r2dxj_x = r2dxj_x, r2yxj_dx = r2yxj_dx, kd = [1, 2, 3], ky = [1, 2, 3], benchmark_covariates = "female")

#

## computes adjusted estimate based on bounds
bound_values = adjusted_estimate(bounds[:, "r2dz_x"], bounds[:, "r2yz_dx"], estimate = 0.0973, se = 0.0232, dof = 783);

## draws contours
ovb_contour_plot(0.097, 0.023, 783);

## adds bounds
add_bound_to_contour(bounds, bound_value = bound_values, label_bump_x = 0.6/30, label_bump_y = 0.6/30)

# ![Figure_9](images/Figure_9.png)
#
# ## References
#
# Cinelli, C. Hazlett, C. (2020) “Making Sense of Sensitivity: Extending Omitted Variable Bias”. Journal of the Royal Statistical Society, Series B (Statistical Methodology). ( [link](https://doi.org/10.1111/rssb.12348) )
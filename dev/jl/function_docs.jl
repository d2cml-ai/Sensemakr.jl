# # Object Documentation
#
# ## `sensemakr` object
# 
# Sensitivity analysis to unobserved confounders.
# 
# The `sensemakr` type performs sensitivity analysis to omitted variables as discussed in Cinelli and Hazlett (2020). It is comprized of several sensitivity statistics for reporting. After creating the object, you may directly use the plot and summary methods of the returned object.
# 
# `sensemakr` is a convenience type. You may use the other sensitivity functions of the package directly, such as the functions for sensitivity plots (`ovb_contour_plot`, `ovb_extreme_plot`), the functions for computing bias-adjusted estimates and t-values (`adjusted_estimate`, `adjusted_t`), the functions for computing the robustness value and partial $R^2$ (`robustness_value`, `partial_r2`), or the functions for bounding the strength of unobserved confounders (`ovb_bounds`), among others.
#
# ```@docs
# Sensemakr.sensemakr
# ```
#
# ```@docs
# summary(::sensemakr, ::Int64)
# ```
#
# ```@docs
# Sensemakr.plot
# ```
#
# ```@docs
# Sensemakr.ovb_minimal_reporting
# ```
# 
# ## Plotting functions
#
# `Sensemakr` provides functions for creating sensitivity contour plots and extreme scenario sensitivity plots. They can be used on an object of class *sensemakr*, directly with an OLS `StatsModels.TableRegressionModel` object, or by providing the required statistics manually.
#
# ```@docs
# Sensemakr.ovb_contour_plot
# ```
#
# ```@docs
# Sensemakr.add_bound_to_contour
# ```
#
# ```@docs
# Sensemakr.ovb_extreme_plot
# ```
#
# ## Sensitivity Bounds
#
# Bounds on the strength of unobserved confounders using observed covariates, as in Cinelli and Hazlett (2020).
#
# These functions may be useful when computing benchmarks for using only summary statistics from papers you see in print.
#
# Currently only the bounds based on partial R2 are implemented. Other bounds will be implemented soon.
#
# ```@docs
# Sensemakr.ovb_bounds
# ```
#
# ```@docs
# Sensemakr.ovb_partial_r2_bound
# ```
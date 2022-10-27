# # Function descriptions
#
# ## Description and analysis functions
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
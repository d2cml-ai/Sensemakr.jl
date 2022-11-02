```@meta
EditURL = "<unknown>/jl/function_docs.jl"
```

# Object Documentation

## `sensemakr` object

Sensitivity analysis to unobserved confounders.

The `sensemakr` type performs sensitivity analysis to omitted variables as discussed in Cinelli and Hazlett (2020). It is comprized of several sensitivity statistics for reporting. After creating the object, you may directly use the plot and summary methods of the returned object.

`sensemakr` is a convenience type. You may use the other sensitivity functions of the package directly, such as the functions for sensitivity plots (`ovb_contour_plot`, `ovb_extreme_plot`), the functions for computing bias-adjusted estimates and t-values (`adjusted_estimate`, `adjusted_t`), the functions for computing the robustness value and partial $R^2$ (`robustness_value`, `partial_r2`), or the functions for bounding the strength of unobserved confounders (`ovb_bounds`), among others.

```@docs
Sensemakr.sensemakr
```

```@docs
summary(::sensemakr, ::Int64)
```

```@docs
Sensemakr.plot
```

```@docs
Sensemakr.ovb_minimal_reporting
```

## Plotting functions

`Sensemakr` provides functions for creating sensitivity contour plots and extreme scenario sensitivity plots. They can be used on an object of class *sensemakr*, directly with an OLS `StatsModels.TableRegressionModel` object, or by providing the required statistics manually.

```@docs
Sensemakr.ovb_contour_plot
```

```@docs
Sensemakr.add_bound_to_contour
```

```@docs
Sensemakr.ovb_extreme_plot
```

## Sensitivity Bounds

Bounds on the strength of unobserved confounders using observed covariates, as in Cinelli and Hazlett (2020).

These functions may be useful when computing benchmarks for using only summary statistics from papers you see in print.

Currently only the bounds based on partial R2 are implemented. Other bounds will be implemented soon.

```@docs
Sensemakr.ovb_bounds
```

```@docs
Sensemakr.ovb_partial_r2_bound
```

## Sensitivity Statistics

Computes the sensitivity statistics: robustness value, partial R2, and Cohen’s f2.

```@docs
Sensemakr.robustness_value
```

```@docs
Sensemakr.partial_r2
```

```@docs
Sensemakr.partial_f2
```

```@docs
Sensemakr.group_partial_r2
```

## Bias functions

Compute bias-adjusted estimates, standard-errors, and t-values.

All methods in the script below have similar purposes and parameters, so they are all described here.

These functions compute bias adjusted estimates (adjusted_estimate), standard-errors (adjusted_se), and t-values (adjusted_t), given a hypothetical strength of the confounder in the partial R2 parameterization.

They return a vector with the adjusted estimate, standard error, or t-value for each partial R2 passed in.

Internally, we also have functions defined to compute the bias and relative_bias, given the same arguments. We also define internal functions to compute the bias function and relative bias function for the partial R2 parameters.

```@docs
Sensemakr.adjusted_estimate
```

```@docs
Sensemakr.adjusted_se
```

```@docs
Sensemakr.adjusted_t
```

```@docs
Sensemakr.adjusted_partial_r2
```

```@docs
Sensemakr.bias
```

````julia
#```@docs
````

Sensemakr.relative_bias
```
## Data

Provides the example data for the package.

```@docs
Sensemakr.load_darfur
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


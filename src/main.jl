"""
    Sensemakr.sensemakr

Object comprised of necessary parameters and statistics for sensitivity analysis
"""
mutable struct sensemakr
    
    model::StatsModels.TableRegressionModel
    treatment::String
    estimate::Float64
    se::Float64
    dof::Int64
    benchmark_covariates::Union{Nothing, String, Array{String}}
    kd::Union{Real, Array{<:Real}}
    ky::Union{Real, Array{<:Real}, Nothing}
    q::Float64
    alpha::Float64
    r2dz_x::Union{Real, Array{<:Real}, Nothing}
    r2yz_dx::Union{Real, Array{<:Real}, Nothing}
    r2dxj_x::Union{Real, Array{<:Real}, Nothing}
    r2yxj_dx::Union{Real, Array{<:Real}, Nothing}
    bound_label
    reduce::Bool
    sensitivity_statistics::Dict
    h0::Float64
    adjusted_estimate::Union{Real, Array{<:Real}, Nothing}
    adjusted_se::Union{Real, Array{<:Real}, Nothing}
    adjusted_t::Union{Real, Array{<:Real}, Nothing}
    adjusted_lower_CI::Union{Real, Array{<:Real}, Nothing}
    adjusted_upper_CI::Union{Real, Array{<:Real}, Nothing}
    bench_bounds::Union{Nothing, DataFrame}
    bounds::Union{Nothing, DataFrame}
end

"""
    Sensemakr.sensemakr(model::StatsModels.TableRegressionModel, treatment::String; kwargs...)

Constructor for the `sensemakr` type. It expects a StatsModels.TableRegressionModel object and a treatment name string as input. It then performs the operations necessary to generate the fields required by the `sensemakr` type.

Arguments:

- `model`: a StatsModels.TableRegressionModel object for the restricted regression model you have provided.
- `treatment`: a string with the name of the "treatment" variable, e.g. the independent variable of interest.
- `benchmark_covariates` (default: `nothing`): a string or vector of strings with the names of the variables to use for benchmark bounding.
- `kd` (default: 1): a float or vector of floats with each being a multiple of the strength of association between a benchmark variable and the treatment variable to test with benchmark bounding.
- `ky` (default: nothing): a float or vector of floats with each being a multiple of the strength of association between a benchmark variable and the treatment variable to test with benchmark bounding.
- `q` (default: 1): a float with the percent to reduce the point estimate by for the robustness value RV_q.
- `alpha` (default = 0.05): a float with the significance level for the robustness value RV_qa to render the estimate not significant.
- `r2dz_x` (default: nothing): a float or vector of floats with the partial R^2 of a putative unobserved confounder "z" with the treatment variable "d", with observed covariates "x" partialed out. In this case, you are manually specifying a putative confounder's strength rather than benchmarking.
- `r2yz_dx` (default: nothing): a float or vector of floats with the partial R^2 of a putative unobserved confounder "z" with the outcome variable "y", with observed covariates "x" and treatment variable "d" partialed out. In this case, you are manually specifying a putative confounder's strength rather than benchmarking.
- `r2dxj_x` (default: nothing): float with the partial R2 of covariate Xj with the treatment D (after partialling out the effect of the remaining covariates X, excluding Xj).
- `r2yxj_dx` (default: nothing): float with the partial R2 of covariate Xj with the outcome Y (after partialling out the effect of the remaining covariates X, excluding Xj).
- `bound_label` (default: "Manual Bound"): a string that specifies what to call the name of a bounding variable, for printing and plotting purposes.
- `reduce` (default: true): whether to reduce (true) or increase (false) the estimate due to putative confounding.
"""
function sensemakr(model::StatsModels.TableRegressionModel, treatment::String; benchmark_covariates = nothing, kd = 1, ky = nothing, q = 1, alpha = 0.05, 
    r2dz_x = nothing, r2yz_dx = nothing, r2dxj_x = nothing, r2yxj_dx = nothing, bound_label = "Manual Bound", reduce = true)

    adjusted_estimate = nothing
    adjusted_se = nothing
    adjusted_t = nothing
    adjusted_lower_CI = nothing
    adjusted_upper_CI = nothing

    if isnothing(ky)
        ky = kd
    end
    if isnothing(r2yz_dx)
        r2yz_dx = r2dz_x
    end
    if isnothing(r2yxj_dx)
        r2yxj_dx = r2dxj_x
    end

    sensitivity_statistics = sensitivity_stats(model = model, treatment = treatment, q = q, alpha = alpha, reduce = reduce)
    estimate = sensitivity_statistics["estimate"][1]
    se = sensitivity_statistics["se"][1]
    dof = Int(dof_residual(model))

    if reduce
        h0 = estimate * (1 - q)
    else
        h0 = estimate * (1 + q)
    end

    if isnothing(r2dz_x)
        bounds = nothing
    else
        r2dz_x, r2yz_dx = check_r2(r2dz_x, r2yz_dx)
        adjusted_estimate = adjusted_estimate(r2dz_x, r2yz_dx, model = model, treatment = treatment, reduce = reduce)
        adjusted_se = adjusted_se(r2dz_x, r2yz_dx; model = model, treatment = treatment)
        adjusted_t = adjusted_t(r2dz_x, r2yz_dx; model = model, treatment = treatment, reduce = reduce)

        se_multiple = quantile(TDist(dof), alpha / 2)
        adjusted_lower_CI = adjusted_estimate - se_multiple * adjusted_se
        adjusted_upper_CI = adjusted_estimate + se_multiple * adjusted_se

        bounds = DataFrame(Dict(
        "r2dz_x" => r2dz_x, 
        "r2yz_dx" => r2yz_dx,
        "bound_label" => repeat([bound_label], length(r2dz_x)), 
        "treatment" => repeat(["treatment"], length(r2dz_x)), 
        "adjusted_estimate" => adjusted_estimate,
        "adjusted_se" => adjusted_se,
        "adjusted_t" => adjusted_t,
        "adjusted_lower_CI" => adjusted_lower_CI,
        "adjusted_upper_CI" => adjusted_upper_CI
        ))
    end

    if !isnothing(benchmark_covariates)
        bench_bounds = ovb_bounds(model, treatment; benchmark_covariates, kd, ky, alpha, h0, reduce)
    elseif !isnothing(r2dxj_x)
        if isnothing(benchmark_covariates)
            benchmark_covariates = "manual_benchmark"
        end
        bench_bounds = obv_partial_r2_bound(r2dxj_x = r2dxj_x, r2yxj_dx = r2yxj_dx, kd = kd, ky = ky, benchmark_covariates = benchmark_covariates)
        if !isnothing(bench_bounds)
            r2dz_x = bench_bounds[:, "r2dz_x"]
            r2yz_dx = bench_bounds[:, "r2yz_dx"]
            bench_bounds[!, "adjusted_estimate"] = adjusted_estimate(r2dz_x, r2yz_dx, estimate = estimate, se = se, dof = dof, reduce = reduce)
            bench_bounds[!, "adjusted_se"] = adjusted_se(r2dz_x, r2yz_dx, se = se, dof = dof)
            bench_bounds[!, "adjusted_estimate"] = adjusted_t(r2dz_x, r2yz_dx, estimate = estimate, se = se, dof = dof, reduce = reduce)
            se_multiple = abs(quantile(TDist(dof), alpha / 2))
            bench_bounds[!, "adjusted_lower_CI"] = bench_bounds[!, "adjusted_estimate"] - se_multiple * bench_bounds[!, "adjusted_se"]
            bench_bounds[!, "adjusted_upper_CI"] = bench_bounds[!, "adjusted_estimate"] + se_multiple * bench_bounds[!, "adjusted_se"]
        end
    else
        bench_bounds = nothing
    end

    if isnothing(bounds)
        bounds = bench_bounds
    else
        append!(bounds, bench_bounds, cols = :union)
    end

    sense_obj = sensemakr(
        model, treatment, estimate, se, dof, benchmark_covariates, kd, ky, q, alpha, r2dz_x, r2yz_dx, r2dxj_x, r2yxj_dx, bound_label, reduce, sensitivity_statistics, 
        h0, adjusted_estimate, adjusted_se, adjusted_t, adjusted_lower_CI, adjusted_upper_CI, bench_bounds, bounds
    )
    return sense_obj
end

#= function sensemakr(estimate::Float64, se::Float64, dof::Int64; benchmark_covariates = nothing, kd = 1, ky = nothing, q = 1, alpha = 0.05, r2dz_x = nothing, 
    r2yz_dx = nothing, r2dxj_x = nothing, r2yxj_dx = nothing, bound_label = "Manual Bound", reduce = true)

    
end =#

"""
    Base.summary(sense_obj::sensemakr, digits::Int64 = 3)

Print a summary of the sensitivity results for `sense_obj`, including robustness value, extreme confounding scenario, and benchmark bounding.

Arguments:

- sense_obj: the sensemakr object to be summarized.
- digits (default: 3): an integer for the number of digits to round numbers to.
- kwargs...: Optional arguments to be dispatched into ovb_contour_plot and ovb_extreme_plot.
"""
function Base.summary(sense_obj::sensemakr, digits::Int64 = 3)

    if sense_obj.reduce
        h0 = round(sense_obj.estimate * (1 - sense_obj.q), digits = digits)
        direction = "reduce"
    else
        h0 = round(sense_obj.estimate * (1 + sense_obj.q), digits = digits)
        direction = "increase"
    end

    println("Sensitivity Analysis to Unobserved Confounding\n")
    if !isnothing(sense_obj.model)
        println("Model Formula: ", sense_obj.model.mf.f, "\n")
    end
    println("Null hypothesis: q = ", sense_obj.q, " and reduce = ", sense_obj.reduce)
    println("-- This means we are considering biases that ", direction, " the absolute value of the current estimate")
    println("-- The null hypothesis deemed problematic is H0:tau = ", h0, "\n")

    println("Unadjusted Estimates of \"", sense_obj.treatment, "\":")
    println("   Coef. Estimate: ", round(sense_obj.estimate, digits = digits))
    println("   Standard Error: ", round(sense_obj.se, digits = digits))
    println("   t-value: ", round(sense_obj.sensitivity_statistics["t_statistic"][1], digits = digits), "\n")

    println("Sensitivity Statistics:")
    println("   Partial R2 of treatment with outcome: ", round.(sense_obj.sensitivity_statistics["r2yd_x"][1], digits = digits))
    println("   Robustness Value, q = ", sense_obj.q, ": ", round(sense_obj.sensitivity_statistics["rv_q"][1], digits = digits))
    println("   Robustness Value, q = ", sense_obj.q, " alpha = ", sense_obj.alpha, ": ", round(sense_obj.sensitivity_statistics["rv_qa"][1], digits = digits), "\n")

    println("Verbal interpretation of sensitivity statistics:\n")
    println(
        "-- Partial R2 of the treatment with the outcome: an extreme confounder (orthogonal to the covariates) ", 
        "that explains 100% of the residual variance of the outcome, would need to explain at least ", 
        round(100.0 * sense_obj.sensitivity_statistics["r2yd_x"][1], digits = digits), " % of the residual variance of the treatment ", 
        "to fully account for the observed estimated effect.\n"
    )

    println(
        "-- Robustness Value, ", "q = ", sense_obj.q, ": unobserved confounders (orthogonal to the covariates) that ",
        "of both the treatment and the outcome are strong enough to bring the point estimate to ", h0,
        " (a bias of ", 100.0 * sense_obj.q, "% of the original estimate). Conversely, unobserved confounders that ",
        "do not explain more than ", round(100.0 * sense_obj.sensitivity_statistics["rv_q"][1], digits = digits), "% of the residual variance ",
        "of both the treatment and the outcome are not strong enough to bring the point estimate to ", h0, ".\n"
    )

    println(
        "-- Robustness Value,", "q = ", sense_obj.q, ", ", "alpha = ", sense_obj.alpha, ": unobserved confounders (orthogonal ",
        "to the covariates) that explain more than ", round(100.0 * sense_obj.sensitivity_statistics["rv_qa"][1], digits = digits), "% of the residual ",
        "variance of both the treatment and the outcome are strong enough to bring the estimate to a range where ",
        "it is no longer 'statistically different' from ", h0, " (a bias of ", 100.0 * sense_obj.q, "% of the original ",
        "estimate), at the significance level of alpha = ", sense_obj.alpha, ". ", "Conversely, unobserved confounders ",
        "that do not explain more than", round(100.0 * sense_obj.sensitivity_statistics["rv_qa"][1], digits = digits), "% of the residual variance",
        "of both the treatment and the outcome are not strong enough to bring the estimate to a range where ",
        "it is no longer 'statistically different' from ", h0, ", at the significance level of alpha = ", sense_obj.alpha, ".\n"
    )

    if !isnothing(sense_obj.bounds)
        println(
            "Bounds on omitted variable bias:\n--The table below shows the maximum strength of unobserved confounders",
            " with association with the treatment and the outcome bounded by a multiple of the observed explanatory",
            " power of the chosen benchmark covariate(s).\n"
        )
        println(sense_obj.bounds)
    end
end

"""
    Sensemakr.plot(sense_obj::sensemakr; kwargs...)

Provide the contour and extreme scenario sensitivity plots of the sensitivity analysis results obtained with the function sensemakr

They are basically dispatchers to the core plot functions ovb_contour_plot and ovb_extreme_plot

Arguments:

- sense_obj: sensemakr obj with the statistics to be plotted.
- plot_type (default: "contour"): Either "extreme" or contour.
- sensitivity_of (default: "estimate"): Either 
"""
function plot(sense_obj::sensemakr; plot_type = "contour", kwargs...)

    if plot_type == "contour"
        ovb_contour_plot(sense_obj; kwargs...)
    elseif plot_type == "extreme"
        ovb_extreme_plot(sense_obj; kwargs...)
    else
        throw(ArgumentError("\"plot_type\" argument must be \"extreme\" or \"contour\""))
    end

end

function ovb_bounds(sense_obj::sensemakr)
    
    model = sense_obj.model
    treatment = sense_obj.treatment
    benchmark_covariates = sense_obj.benchmark_covariates
    kd = sense_obj.kd
    ky = sense_obj.ky
    alpha = sense_obj.alpha
    h0 = sense_obj.h0
    reduce = sense_obj.reduce
    bound = sense_obj.bound
    adjusted_estimates = sense_obj.adjusted_estimates

    return ovb_bounds(model, treatment; benchmark_covariates, kd, ky, alpha, h0, reduce, bound, adjusted_estimates)
end

"""
    Base.print(sense_obj::sensemakr, digits = 3)

Print a short summary of the sensitivity results for a sensemakr object, including formula, hypothesis, and sensitivity analysis.

Arguments:

- `sense_obj`: the sensemakr object to be summarized.
- `digits` (default: 3): an integer for the number of digits to round numbers to.
"""
function Base.print(sense_obj::sensemakr, digits = 3)

    if sense_obj.reduce
        h0 = round(sense_obj.estimate * (1 - sense_obj.q), digits = digits)
        direction = "reduce"
    else
        h0 = round(sense_obj.estimate * (1 + sense_obj.q), digits = digits)
        direction = "increase"
    end

    println("Sensitivity Analysis to Unobserved Confounding\n")
    if !isnothing(sense_obj.model)
        println("Model Formula: ", sense_obj.model.mf.f, "\n")
    end
    println("Null hypothesis: q = ", sense_obj.q, " and reduce = ", sense_obj.reduce, "\n")

    println("Unadjusted Estimates of \"", sense_obj.treatment, "\":")
    println("   Coef. Estimate: ", round(sense_obj.estimate, digits = digits))
    println("   Standard Error: ", round(sense_obj.se, digits = digits))
    println("   t-value: ", round(sense_obj.sensitivity_statistics["t_statistic"][1], digits = digits), "\n")

    println("Sensitivity Statistics:")
    println("   Partial R2 of treatment with outcome: ", round.(sense_obj.sensitivity_statistics["r2yd_x"][1], digits = digits))
    println("   Robustness Value, q = ", sense_obj.q, ": ", round(sense_obj.sensitivity_statistics["rv_q"][1], digits = digits))
    println("   Robustness Value, q = ", sense_obj.q, " alpha = ", sense_obj.alpha, ": ", round(sense_obj.sensitivity_statistics["rv_qa"][1], digits = digits), "\n")
end

"""
    ovb_minimal_reporting(sense_obj::sensemakr, digits::Int64 = 3, res_display::Bool = true)

This function returns the HTML code for a table summarizing the sensemakr object.

Arguments:

- `digits` (default: 3): an integer for the number of digits to round numbers to.
- `res_display` (default: true): whether to display the table.
"""
function ovb_minimal_reporting(sense_obj::sensemakr, digits::Int64 = 3, res_display::Bool = true)

    result = "<table style = 'align:center'>\n" * "<thead\n" * 
    "<tr>" * 
    "\t<th style='text-align:center;border-bottom: 1px solid black;border-top: 1px solid black> </th>\n" * 
    "\t<th colspan = 6 style='text-align:center;border-bottom: 1px solid black; border-top: 1px solid black'> Outcome: " * 
    string(sense_obj.model.mf.f.lhs) * "</tr>\n" * 
    "</tr>\n" * 
    "<tr>\n" * 
    "\t<th style='text-align:left;border-top: 1px solid black'> Treatment </th>\n" * 
    "\t<th style='text-align:right;border-top: 1px solid black'> Est. </th>\n" * 
    "\t<th style='text-align:right;border-top: 1px solid black'> S.E. </th>\n" * 
    "\t<th style='text-align:right;border-top: 1px solid black'> t-value </th>\n" * 
    "\t<th style='text-align:right;border-top: 1px solid black'> R<sup>2</sup><sub>Y~D|X</sub> </th>\n" * 
    "\t<th style='text-align:right;border-top: 1px solid black'> RV<sub>q = " * 
    string(sense_obj.q) * "</sub> </th>\n" * 
    "\t<th style='text-align:right;border-top: 1px solid black'> RV<sub>q = " * 
    string(sense_obj.q) * ", &alpha; = " * 
    string(sense_obj.alpha) * "</sub> </th>\n" * 
    "</tr>\n" * 
    "</thead>\n" * 
    "<tbody>\n <tr> \n" * 
    "\t<td style 'text-alighn:left;border-bottom: 1px solid black'><i>" * 
    string(sense_obj.treatment) * "</i></td>\n" * 
    "\t<td style='text-align:right;border-bottom: 1px solid black'>" * 
    string(round(sense_obj.sensitivity_statistics["estimate"][1], digits = digits)) * " </td>\n" *
    "\t<td style='text-align:right;border-bottom: 1px solid black'>" * 
    string(round(sense_obj.sensitivity_statistics["se"][1], digits = digits)) * " </td>\n" * 
    "\t<td style='text-align:right;border-bottom: 1px solid black'>" * 
    string(round(sense_obj.sensitivity_statistics["t_statistic"][1], digits = digits)) * " </td>\n" * 
    "\t<td style='text-align:right;border-bottom: 1px solid black'>" * 
    string(round(sense_obj.sensitivity_statistics["r2yd_x"][1] * 100, digits = digits - 2)) * "% </td>\n" * 
    "\t<td style='text-align:right;border-bottom: 1px solid black'>" * 
    string(round(sense_obj.sensitivity_statistics["rv_q"][1] * 100, digits = digits - 2)) * "% </td>\n" * 
    "\t<td style='text-align:right;border-bottom: 1px solid black'>" * 
    string(round(sense_obj.sensitivity_statistics["rv_qa"][1] * 100, digits = digits - 2)) * "% </td>\n" *
    "</tr>\n</tbody>\n" 
    
    if isnothing(sense_obj.bounds)
        result = result * 
        "<tr>\n" * 
        "<td colspan = 7 style='text-align:right;border-top: 1px solid black;border-bottom: 1px solid transparent;font-size:11px'>" * 
        "Note: df = " * string(sense_obj.sensitivity_statistics["dof"]) *  
        "</td>\n" * 
        "</tr>\n" * 
        "</table>"
    else
        result = result * "<tr>\n" * 
        "<td colspan = 7 style='text-align:right;border-top: 1px solid black;border-bottom: 1px solid transparent;font-size:11px'>" * 
        "Note: df = " * string(sense_obj.sensitivity_statistics["dof"]) * "; " * 
        "Bound ( " * string(sense_obj.bounds[:, "bound_label"][1]) * " ):  " * 
        "R<sup>2</sup><sub>Y~Z|X,D</sub> =  " * 
        string(round(sense_obj.bounds[:, "r2yz_dx"][1] * 100, digits = digits - 2)) * 
        "%, R<sup>2</sup><sub>D~Z|X</sub> =" * 
        string(round(sense_obj.bounds[:, "r2dz_x"][1] * 100, digits = digits - 2)) * 
        "%" *  
        "</td>\n" * 
        "</tr>\n" * 
        "</table>"
    end

    
    if res_display
        result = HTML(result)
        return result
    else
        return result
    end
end

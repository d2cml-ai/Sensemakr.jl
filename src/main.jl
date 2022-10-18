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
        bench_bounds = obv_bounds(model, treatment; benchmark_covariates, kd, ky, alpha, h0, reduce)
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
        println("Model Formula: ", string.(keys(sense_obj.model.mf.schema))[1], " ~ ", join(string.(keys(sense_obj.model.mf.schema))[2:end], " + "), "\n")
    end
    println("Null hypothesis: q = ", sense_obj.q, " and reduce = ", sense_obj.reduce)
    println("-- This means we are considering biases that ", direction, " the absolute value of the current estimate")
    println("-- The null hypothesis deemed problematic is H0:tau = ", h0, "\n")

    println("Unadjusted Estimates of \"", sense_obj.treatment, "\":")
    println("   Coef. Estimate: ", round(sense_obj.estimate, digits = digits))
    println("   Standard Error: ", round(sense_obj.se, digits = digits))
    println("   t-value: ", round(sense_obj.sensitivity_statistics["t_statistic"][1], digits = digits))

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

function plot(sense_obj::sensemakr; plot_type = "contour", kwargs...)

    if plot_type == "contour"
        ovb_contour_plot(sense_obj; kwargs...)
    elseif plot_type == "extreme"
        ovb_extreme_plot(sense_obj; kwargs...)
    else
        throw(ArgumentError("\"plot_type\" argument must be \"extreme\" or \"contour\""))
    end

end

function obv_bounds(sense_obj::sensemakr)
    
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

    return obv_bounds(model, treatment; benchmark_covariates, kd, ky, alpha, h0, reduce, bound, adjusted_estimates)
end

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
        println("Model Formula: ", string.(keys(sense_obj.model.mf.schema))[1], " ~ ", join(string.(keys(sense_obj.model.mf.schema))[2:end], " + "), "\n")
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


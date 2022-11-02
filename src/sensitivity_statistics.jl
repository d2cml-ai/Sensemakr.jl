"""
    robustness_value(; model::Union{Nothing, StatsModels.TableRegressionModel} = nothing, covariates::Union{Nothing, String, Vector{String}} = nothing, t_statistic::Union{Real, Nothing} = nothing, dof::Union{Nothing, Int64} = nothing, kwargs...)

Compute the robustness value of a regression coefficient.

The robustness value describes the minimum strength of association (parameterized in terms of partial R2) that omitted variables would need to have both with the treatment and with the outcome to change the estimated coefficient by a certain amount (for instance, to bring it down to zero).

For instance, a robustness value of 1% means that an unobserved confounder that explain 1% of the residual variance of the outcome and 1% of the residual variance of the treatment is strong enough to explain away the estimated effect. Whereas a robustness value of 90% means that any unobserved confounder that explain less than 90% of the residual variance of both the outcome and the treatment assignment cannot fully account for the observed effect. You may also compute robustness value taking into account sampling uncertainty. See details in Cinelli and Hazlett (2020).

# Arguments:

- `model`: object for the restricted regression you have provided.
- `covariates`: names of the variables to use for benchmark bounding.
- `t_statistic`: the t_statistic for the restricted model regression.
- `dof`: degrees of freedom of the restricted regression
- `q` (default: 1): a float with the percent to reduce the point estimate by for the robustness value RV_q.
- `alpha` (default = 1.0): a float with the significance level for the robustness value RV_qa to render the estimate not significant.
"""
function robustness_value(; model::Union{Nothing, StatsModels.TableRegressionModel} = nothing, covariates::Union{Nothing, String, Vector{String}} = nothing, 
    t_statistic::Union{Real, Nothing} = nothing, dof::Union{Nothing, Int64} = nothing, q::Real = 1, alpha::Float64 = 1.0)
    
    if isnothing(model) && (isnothing(t_statistic) || isnothing(dof))
        throw(ArgumentError("robustness_value requires either a GLM LinearModel object or a t-statistic and degrees of freedom"))
    end
    check_q(q)
    check_alpha(alpha)
    
    if !isnothing(model)
        model_data = model_helper(model, covariates)
        t_statistic = model_data["t_statistics"]
        dof = model_data["dof"]
    elseif t_statistic isa Real
        t_statistic = [t_statistic]
    end

    fq = q * abs.(t_statistic / sqrt(dof))
    f_crit = abs(quantile(TDist(dof - 1), alpha / 2)) / sqrt(dof - 1)
    fqa = fq .- f_crit

    rv = 0.5 * (sqrt.(fqa .^ 4 + (4 * fqa .^ 2)) .- fqa .^ 2)
    rvx = (fq .^ 2 .- f_crit ^ 2) ./ (1 .+ fq .^ 2)

    rv_out = rv
    
    rv_out[fqa .< 0] .= 0
    rv_out[((fqa .> 0) .&& (fq .> (1 / f_crit)))] = rvx[((fqa .> 0) .&& (fq .> (1 / f_crit)))]
    
    return rv_out
end

"""
    partial_r2(; model::Union{Nothing, StatsModels.TableRegressionModel} = nothing, covariates::Union{Nothing, String, Vector{String}} = nothing, t_statistic::Union{Real, Nothing} = nothing, dof::Union{Nothing, Int64} = nothing)

Compute the partial R2 for a linear regression model.

The partial R2 describes how much of the residual variance of the outcome (after partialing out the other covariates) a covariate explains.

The partial R2 can be used as an extreme-scenario sensitivity analysis to omitted variables. Considering an unobserved confounder that explains 100% of the residual variance of the outcome, the partial R2 describes how strongly associated with the treatment this unobserved confounder would need to be in order to explain away the estimated effect.

For details see Cinelli and Hazlett (2020).

# Arguments:

- `model`: object for the restricted regression you have provided.
- `covariates`: names of the variables to use for benchmark bounding.
- `t_statistic`: the t_statistic for the restricted model regression.
- `dof`: degrees of freedom of the restricted regression
"""
function partial_r2(; model::Union{Nothing, StatsModels.TableRegressionModel} = nothing, covariates::Union{Nothing, String, Vector{String}} = nothing, 
    t_statistic::Union{Vector{<:Real}, Real, Nothing} = nothing, dof::Union{Nothing, Int64} = nothing)
    
    if isnothing(model) && (isnothing(t_statistic) || isnothing(dof))
        throw(ArgumentError("partial_r2 requires either a GLM LinearModel object or a t-statistic and degrees of freedom"))
    end
    
    if !isnothing(model)
        model_data = model_helper(model, covariates)
        t_statistic = model_data["t_statistics"]
        dof = model_data["dof"]
    end

    return t_statistic .^ 2 ./ (t_statistic .^ 2 .+ dof)
end

"""
    partial_f2(; model::Union{Nothing, StatsModels.TableRegressionModel} = nothing, covariates::Union{Nothing, String, Vector{String}} = nothing, t_statistic::Union{Real, Nothing} = nothing, dof::Union{Nothing, Int64} = nothing)

Compute the partial (Cohen's) f2 for a linear regression model.

The partial (Cohen's) f2 is a common measure of effect size (a transformation of the partial R2) that can also be used directly for sensitivity analysis using a bias factor table. For details see Cinelli and Hazlett (2020).

# Arguments

- `model`: object for the restricted regression you have provided.
- `covariates`: names of the variables to use for benchmark bounding.
- `t_statistic`: the t_statistic for the restricted model regression.
- `dof`: degrees of freedom of the restricted regression
"""
function partial_f2(; model::Union{Nothing, StatsModels.TableRegressionModel} = nothing, covariates::Union{Nothing, String, Vector{String}} = nothing, 
    t_statistic::Union{Vector{<:Real}, Real, Nothing} = nothing, dof::Union{Nothing, Int64} = nothing)

    if isnothing(model) && (isnothing(t_statistic) || isnothing(dof))
        throw(ArgumentError("partial_f2 requires either a GLM LinearModel object or a t-statistic and degrees of freedom"))
    end

    if !isnothing(model)
        model_data = model_helper(model, covariates)
        t_statistic = model_data["t_statistics"]
        dof = model_data["dof"]
    end

    return t_statistic .^ 2 / dof
end

function partial_f(; model = nothing, covariates = nothing, t_statistic = nothing, dof = nothing)

    return sqrt.(partial_f2(model = model, covariates = covariates, t_statistic = t_statistic, dof = dof))
end

"""
Partial R2 of groups of covariates in a linear regression model.

This function computes the partial R2 of a group of covariates in a linear regression model. Multivariate version of the partial_r2 function; see that for more details.

# Arguments

- `model`: object for the restricted regression you have provided.
- `covariates`: names of the variables to use for benchmark bounding.
- `t_statistic`: t_statistic for the restricted model regression.
- `p`: number of parameters.
- `dof`: degrees of freedom of the restricted regression.
"""
function group_partial_r2(; model::Union{Nothing, StatsModels.TableRegressionModel} = nothing, covariates::Union{Nothing, String, Vector{String}} = nothing, f_statistic::Union{Vector{<:Real}, Real, Nothing} = nothing, p::Union{Int64, Nothing} = nothing, dof::Union{Int64, Nothing} = nothing)

    if (isnothing(model) || isnothing(covariates)) && (isnothing(f_statistic) || isnothing(p) || isnothing(dof))
        throw(ArgumentError("group_partial_r2 requires either a GLM LinearModel object and covariates or an f-statistic, number of parameters, and degrees of freedom"))
    end

    if isnothing(f_statistic) || isnothing(p) || isnothing(dof)
        params = coef(model)
        check_covariates(coefnames(model), covariates)
        var_index = findall(in.(coefnames(model), Ref(covariates)))
        params = params[var_index]
        if size(params, 1) == 1
            return partial_r2(model = model, covariates = covariates, t_statistic = f_statistic, dof = dof)
        end
        v = vcov(model)[var_index, var_index]
        dof = dof_residual(model)
        p = size(params, 1)
        f_statistic = (params' * inv(v)) * params / p
    end
    return [f_statistic * p / (f_statistic * p + dof)]
end

"""
    sensitivity_stats(; model::Union{Nothing, StatsModels.TableRegressionModel} = nothing, treatment::Union{String, Nothing} = nothing, estimate::Union{Nothing, Real} = nothing, se::Union{Real, Nothing} = nothing, dof::Union{Real, Nothing} = nothing)

Computes the robustness_value, partial_r2 and partial_f2 of the coefficient of interest.

# Arguments:

- `model`: object for the restricted regression you have provided.
- `treatment`: names of the treatment variable.
- `estimate`: coefficient estimate of the restricted regression.
- `se`: standard error of the restricted regression.
- `dof`: degrees of freedom of the restricted regression.
- `q` (default: 1): a float with the percent to reduce the point estimate by for the robustness value RV_q.
- `alpha` (default = 1.0): a float with the significance level for the robustness value RV_qa to render the estimate not significant.
- `reduce` (default = true): whether to reduce or increase the estimate due to confounding.
"""
function sensitivity_stats(; model::Union{Nothing, StatsModels.TableRegressionModel} = nothing, treatment::Union{String, Nothing} = nothing, estimate::Union{Nothing, Real} = nothing, se::Union{Real, Nothing} = nothing, 
    dof::Union{Real, Nothing} = nothing, q::Real = 1, alpha::Float64 = 0.05, reduce::Bool = true)

    if (isnothing(model) || isnothing(treatment)) && (isnothing(estimate) || isnothing(se) || isnothing(dof))
        throw(ArgumentError("sensitivity_stats requires either a GLM LinearModel object and treatment name or an f-statistic, number of parameters, and degrees of freedom"))
    end
    if !isnothing(model)
        model_data = model_helper(model, treatment)
        estimate = model_data["estimate"]
        se = model_data["se"]
        dof = Int(model_data["dof"])
    end

    check_q(q)
    check_alpha(alpha)
    check_se(se)
    check_dof(dof)
    
    if reduce
        h0 = estimate * (1 - q)
    else
        h0 = estimate * (1 + q)
    end
    original_t = estimate ./ se
    t_statistic = (estimate - h0) ./ se
    r2yd_x = partial_r2(t_statistic = original_t, dof = dof)
    rv_q = robustness_value(t_statistic = original_t, dof = dof, q = q)
    rv_qa = robustness_value(t_statistic = original_t, dof = dof, q = q, alpha = alpha)
    f2yd_x = partial_f2(t_statistic = original_t, dof = dof)
    sensitivity_stats_df = Dict(
        "estimate" => estimate,
        "se" => se,
        "t_statistic" => t_statistic,
        "r2yd_x" => r2yd_x, 
        "rv_q" => rv_q,
        "rv_qa" => rv_qa,
        "f2yd_x" => f2yd_x,
        "dof" => dof
    )
    return sensitivity_stats_df
end

function model_helper(model, covariates = nothing)
    
    error_if_no_dof(model)
    if !isnothing(covariates)
        covariates = check_covariates(coefnames(model), covariates)
        used_variables = covariates
    else
        used_variables = coefnames(model)
    end
    var_index = findall(in.(coefnames(model), Ref(used_variables)))
    model_info = Dict(
        "covariates" => used_variables,
        "estimate" => coef(model)[var_index],
        "se" => stderror(model)[var_index],
        "t_statistics" => coeftable(model).cols[3][var_index],
        "dof" => Int(dof_residual(model))
    )
    return model_info
end

function check_r2(r2dz_x::Union{Nothing, Real, Array{Int64}, Array{Float64}}, r2yz_dx::Union{Real, Array{Int64}, Array{Float64}})

    if isnothing(r2dz_x)
        return r2dz_x, r2yz_dx
    end
    for r2s in [r2dz_x, r2yz_dx]
        if count(i -> !(0 <= i <= 1), r2s) > 0
            throw(DomainError(i, "partial R2 must be number or array of numbers in [0, 1]"))
        end
    end
    return r2dz_x, r2yz_dx
end

function check_q(q::Real)
    
    if q < 0
        throw(DomainError(q, "q must be greater than 0"))
    end
end

function check_alpha(alpha::Float64)

    if !(0.0 <= alpha <= 1.0)
        throw(DomainError(alpha, "alpha must be in [0, 1]"))
    end
end

function check_se(se::Union{Real, Vector})
    
    if se isa Vector
        if length(se) > 1
            throw(ArgumentError("se must contain a single number"))
        end
        se = se[1]
    end
    if se < 0
        throw(DomainError(se, "se must be greater than 0"))
    end
end

function check_dof(dof::Real)
    
    if dof <= 0
        throw(DomainError(dof, "degrees of freedom must be greater than 0"))
    end
    dof = Int(dof)
end

function error_if_no_dof(model)
    
    if dof_residual(model) == 0
        throw(ErrorException("There are 0 residual degrees of freedom in the regression provided"))
    end
end

function check_covariates(all_names, covariates)
    if !isnothing(covariates)
        if covariates isa String
            covariates = [covariates]
        end
        if 0 in isa.(covariates, String)
            throw(TypeError)
        end
        not_found = covariates[(!in).(covariates, Ref(all_names))]
        if size(not_found, 1) > 0
            throw(ErrorException("Variables not found in model: " * join(not_found, ", ")))
        end
    end
    return covariates
end
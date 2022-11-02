"""
    adjusted_estimate(r2dz_x::Union{Float64, Vector{Float64}}, r2yz_dx::Union{Float64, Vector{Float64}}; kwargs...)

Compute the bias-adjusted coefficient estimate.

# Arguments

- `r2dz_x`: partial \$R^2\$ of a putative unobserved confounder "z" with the treatment variable "d", with observed covariates "x" partialed out.
- `r2yz_dx`: partial \$R^2\$ of a putative unobserved confounder "z" with the outcome "y", with observed covariates "x" and treatment "d" partialed out.
- `model`: `StatsModels.TableRegressionModel` object for the restricted regression you have provided.
- `treatment`: string with the name of the treatment variable of interest.
- `estimate`: float with the unadjusted estimate of the coefficient for the independent variable of interest.
- `se`: float with the unadjusted standard error of the regression.
- `dof`: an int with the degrees of freedom of the regression.
- `reduce` (default: true): whether to reduce (true) or increase (false) the estimate due to putative confounding.
"""
function adjusted_estimate(r2dz_x::Union{Float64, Vector{Float64}}, r2yz_dx::Union{Float64, Vector{Float64}}; model = nothing, treatment = nothing, estimate = nothing, se = nothing, 
    dof = nothing, reduce::Bool = true)

    r2dz_x, r2yz_dx, estimate, se, dof = param_check("adjusted_estimate", r2dz_x, r2yz_dx, model = model, treatment = treatment, estimate = estimate, se = se, dof = dof)

    if reduce
        return sign.(estimate) .* (abs.(estimate) .- bias(r2dz_x, r2yz_dx, se = se, dof = dof))
    else
        return sign.(estimate) .* (abs.(estimate) .+ bias(r2dz_x, r2yz_dx, se = se, dof = dof))
    end
end

"""
    adjusted_se(r2dz_x::Union{Float64, Vector{Float64}}, r2yz_dx::Union{Float64, Vector{Float64}}; kwargs...)

Compute the bias-adjusted regression standard error.

# Arguments

- `r2dz_x`: partial \$R^2\$ of a putative unobserved confounder "z" with the treatment variable "d", with observed covariates "x" partialed out.
- `r2yz_dx`: partial \$R^2\$ of a putative unobserved confounder "z" with the outcome "y", with observed covariates "x" and treatment "d" partialed out.
- `model`: `StatsModels.TableRegressionModel` object for the restricted regression you have provided.
- `treatment`: string with the name of the treatment variable of interest.
- `se`: float with the unadjusted standard error of the regression.
- `dof`: an int with the degrees of freedom of the regression.
"""
function adjusted_se(r2dz_x::Union{Float64, Vector{Float64}}, r2yz_dx::Union{Float64, Vector{Float64}}; model = nothing, treatment = nothing, se = nothing, dof = nothing)

    r2dz_x, r2yz_dx, estimate, se, dof = param_check("adjusted_se", r2dz_x, r2yz_dx, model = model, treatment = treatment, se = se, dof = dof, estimate_is_param = false)
    new_se = sqrt.((1 .- r2yz_dx) ./ (1 .- r2dz_x)) .* se .* sqrt(dof / (dof - 1))
    return new_se
end

"""
    adjusted_t(r2dz_x::Union{Float64, Vector{Float64}}, r2yz_dx::Union{Float64, Vector{Float64}}; kwargs...)

Compute bias-adjusted t-statistic, (adjusted_estimate - h0) / adjusted_se.

# Arguments

- `r2dz_x`: partial \$R^2\$ of a putative unobserved confounder "z" with the treatment variable "d", with observed covariates "x" partialed out.
- `r2yz_dx`: partial \$R^2\$ of a putative unobserved confounder "z" with the outcome "y", with observed covariates "x" and treatment "d" partialed out.
- `model`: `StatsModels.TableRegressionModel` object for the restricted regression you have provided.
- `treatment`: string with the name of the treatment variable of interest.
- `estimate`: float with the unadjusted estimate of the coefficient for the independent variable of interest.
- `se`: float with the unadjusted standard error of the regression.
- `dof`: an int with the degrees of freedom of the regression.
- `reduce` (default: true): whether to reduce (true) or increase (false) the estimate due to putative confounding.
- `h0`: (default: 0): the test value for null hypothesis.
"""
function adjusted_t(r2dz_x::Union{Float64, Vector{Float64}}, r2yz_dx::Union{Float64, Vector{Float64}}; model = nothing, treatment = nothing, estimate = nothing, se = nothing, 
    dof = nothing, reduce = true, h0::Real = 0)

    r2dz_x, r2yz_dx, estimate, se, dof = param_check("adjusted_t", r2dz_x, r2yz_dx, model = model, treatment = treatment, estimate = estimate, se = se, dof = dof)
    new_estimate = adjusted_estimate(r2dz_x, r2yz_dx, estimate = estimate[1], se = se, dof = dof, reduce = reduce)
    new_t = (new_estimate .- h0) ./ adjusted_se(r2dz_x, r2yz_dx, se = se, dof = dof)
    return new_t
end


"""
    adjusted_partial_r2(r2dz_x::Union{Float64, Vector{Float64}}, r2yz_dx::Union{Float64, Vector{Float64}}; kwargs...)

Compute the bias-adjusted partial R2, based on adjusted_t

# Arguments

- `r2dz_x`: partial \$R^2\$ of a putative unobserved confounder "z" with the treatment variable "d", with observed covariates "x" partialed out.
- `r2yz_dx`: partial \$R^2\$ of a putative unobserved confounder "z" with the outcome "y", with observed covariates "x" and treatment "d" partialed out.
- `model`: `StatsModels.TableRegressionModel` object for the restricted regression you have provided.
- `treatment`: string with the name of the treatment variable of interest.
- `estimate`: float with the unadjusted estimate of the coefficient for the independent variable of interest.
- `se`: float with the unadjusted standard error of the regression.
- `dof`: an int with the degrees of freedom of the regression.
- `reduce` (default: true): whether to reduce (true) or increase (false) the estimate due to putative confounding.
- `h0`: (default: 0): the test value for null hypothesis.
"""
function adjusted_partial_r2(r2dz_x::Union{Float64, Vector{Float64}}, r2yz_dx::Union{Float64, Vector{Float64}}; model = nothing, treatment = nothing, estimate = nothing, 
    se = nothing, dof = nothing, reduce = true, h0 = 0)

    r2dz_x, r2yz_dx, estimate, se, dof = param_check("adjusted_partial_r2", r2dz_x, r2yz_dx, model = model, treatment = treatment, estimate = estimate, se = se, dof = dof)
    new_t = adjusted_t(r2dz_x, r2yz_dx, estimate = estimate, se = se, dof = dof, reduce = reduce, h0 = h0)
    return partial_r2(t_statistic = new_t, dof = (dof - 1))
end

"""
    bias(r2dz_x::Union{Float64, Vector{Float64}}, r2yz_dx::Union{Float64, Vector{Float64}}; kwargs...)

Compute the omitted variable bias for the partial R2 parameterization.

# Arguments

- `r2dz_x`: partial \$R^2\$ of a putative unobserved confounder "z" with the treatment variable "d", with observed covariates "x" partialed out.
- `r2yz_dx`: partial \$R^2\$ of a putative unobserved confounder "z" with the outcome "y", with observed covariates "x" and treatment "d" partialed out.
- `model`: `StatsModels.TableRegressionModel` object for the restricted regression you have provided.
- `treatment`: string with the name of the treatment variable of interest.
- `se`: float with the unadjusted standard error of the regression.
- `dof`: an int with the degrees of freedom of the regression.
"""
function bias(r2dz_x::Union{Float64, Vector{Float64}}, r2yz_dx::Union{Float64, Vector{Float64}}; model = nothing, treatment = nothing, se = nothing, dof = nothing)

    r2dz_x, r2yz_dx, estimate, se, dof = param_check("bias", r2dz_x, r2yz_dx, model = model, treatment = treatment, se = se, dof = dof, estimate_is_param = false)
    bias_val = bf(r2dz_x, r2yz_dx) .* se .* sqrt(dof)
    return bias_val
end

"""
    relative_bias(r2dz_x::Union{Float64, Vector{Float64}}, r2yz_dx::Union{Float64, Vector{Float64}}; kwargs...)

Compute the relative bias for the partial R2 parameterization.

# Arguments

- `r2dz_x`: partial \$R^2\$ of a putative unobserved confounder "z" with the treatment variable "d", with observed covariates "x" partialed out.
- `r2yz_dx`: partial \$R^2\$ of a putative unobserved confounder "z" with the outcome "y", with observed covariates "x" and treatment "d" partialed out.
- `model`: `StatsModels.TableRegressionModel` object for the restricted regression you have provided.
- `treatment`: string with the name of the treatment variable of interest.
- `estimate`: float with the unadjusted estimate of the coefficient for the independent variable of interest.
- `se`: float with the unadjusted standard error of the regression.
- `dof`: an int with the degrees of freedom of the regression.
"""
function relative_bias(r2dz_x::Union{Float64, Vector{Float64}}, r2yz_dx::Union{Float64, Vector{Float64}}; model = nothing, treatment = nothing, estimate = nothing, se = nothing, 
    dof = nothing)
    
    r2dz_x, r2yz_dx, estimate, se, dof = param_check("relative_bias", r2dz_x, r2yz_dx, model = model, treatment = treatment, estimate = estimate, se = se, dof = dof)
    t_statistic = abs(estimate / se)
    f = partial_f(t_statistic = t_statistic, dof =  dof)
    bf_val = bf(r2dz_x, r2yz_dx)
    q = bf_val ./ f
    return q
end

function rel_bias(r_est, est)

    return (r_est .- est) ./ r_est
end

function bf(r2dz_x, r2yz_dx)

    return sqrt.((r2yz_dx .* r2dz_x) ./ (1 .- r2dz_x))
end

function param_check(function_name::String, r2dz_x, r2yz_dx; model = nothing, treatment::Union{String, Nothing} = nothing, estimate::Union{Nothing, Real} = nothing, se = nothing, 
    dof = nothing, estimate_is_param::Bool = true)

    if estimate_is_param
        if (isnothing(model) || isnothing(treatment)) && (isnothing(estimate) || isnothing(se) || isnothing(dof))
            throw(ArgumentError(join(["in addition to r2dz_x and r2yz_x, ", 
            function_name, 
            " requires either a GLM LinearModel and a treatment name or the current estimate, standard error, and degrees of freedom"])))
        end
    else
        if (isnothing(model) || isnothing(treatment)) && (isnothing(se) || isnothing(dof))
            throw(ArgumentError(join(["in addition to r2dz_x and r2yz_x, ", 
            function_name, 
            " requires either a GLM LinearModel and a treatment name or the current standard error and degrees of freedom"])))
        end
    end

    if !isnothing(model)
        model_data = model_helper(model, treatment)
        estimate = model_data["estimate"]
        se = model_data["se"]
        dof = model_data["dof"]
    end
    check_se(se)
    check_dof(dof)
    r2dz_x, r2yz_dx = check_r2(r2dz_x, r2yz_dx)
    return r2dz_x, r2yz_dx, estimate, se, dof
end
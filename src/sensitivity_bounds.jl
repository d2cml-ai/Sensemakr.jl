"""
    ovb_bounds(model::StatsModels.TableRegressionModel, treatment::String; kwargs...)

Returns a `DataFrame` which provides bounds on the strength of unobserved confounders using observed covariates, as in Cinelli and Hazlett (2020).

# Arguments

- `model`: object for the restricted regression you have provided.
- `treatment`: name of the treatment variable or variable of interest.
- `benchmark_covariates`: names of variables to use for benchmark bounding.
- `kd` (default: 1): multiple of the strength of association between a benchmark variable and the treatment.
- `ky`: same as kd except measured in terms of strength of association with the outcome variable.
- `alpha` (default: 0.05): significance level for the robustness value \$RV_qa\$ to render the estimate not significant.
- `h0` (default: 0): null hypothesis effect size.
- `reduce` (default: true): whether to reduce (true) or increase (false) the estimate due to putative computing.
- `bound` (default: "partial_r2"): type of bound to perform; as of now, only partial R^2 bounding is allowed.
- `adjusted_estimates` (default: true): whether to compute bias-adjusted estimates, standard errors, and t-statistics.
"""
function ovb_bounds(model::StatsModels.TableRegressionModel, treatment::String; benchmark_covariates::Union{String, Vector{String}, Nothing} = nothing, 
    kd::Union{Vector{<:Real}, Real} = 1, ky::Union{Vector{<:Real}, Real, Nothing} = nothing, alpha::Float64 = 0.05, h0::Real = 0, reduce::Bool = true, 
    bound::String = "partial_r2", adjusted_estimates::Bool = true)

    if isnothing(ky)
        ky = kd
    end
    if bound != "partial_r2"
        throw(ErrorException("only partial r2 implemented as of now"))
    end
    if benchmark_covariates isa String
        benchmark_covariates = [benchmark_covariates]
    end
    bounds = ovb_partial_r2_bound(model = model, treatment = treatment, benchmark_covariates = benchmark_covariates, kd = kd, ky = ky)
    
    if adjusted_estimates
        bounds[!, "treatment"] = repeat([treatment], length(kd) * length(benchmark_covariates))
        bounds[!, "adjusted_estimate"] = adjusted_estimate(bounds[:, "r2dz_x"], bounds[:, "r2yz_dx"], model = model, treatment = treatment, reduce = reduce)
        bounds[!, "adjusted_se"] = adjusted_se(bounds[:, "r2dz_x"], bounds[:, "r2yz_dx"], model = model, treatment = treatment)
        bounds[!, "adjusted_t"] = adjusted_t(bounds[:, "r2dz_x"], bounds[:, "r2yz_dx"], model = model, treatment = treatment, reduce = reduce, h0 = h0)

        se_multiple = abs(quantile(TDist(dof_residual(model)), alpha / 2))
        bounds[!, "adjusted_lower_CI"] = bounds[:, "adjusted_estimate"] - se_multiple * bounds[:, "adjusted_se"]
        bounds[!, "adjusted_upper_CI"] = bounds[:, "adjusted_estimate"] + se_multiple * bounds[:, "adjusted_se"]
    end
    return bounds
end

"""
    ovb_partial_r2_bound(; kwargs...)

Returns a `DataFrame` with the bounds on the strength of the unobserved confounder.

Adjusted estimates, standard errors and t-values (among other quantities) need to be computed manually by the user using those bounds with the functions adjusted_estimate, adjusted_se and adjusted_t.

# Arguments

- `model`: `StatsModels.TableRegressionModel` object for the restricted regression you have provided.
- `treatment`: string with the name of the treatment variable of interest.
- `r2dxj_x`: float with the partial R2 of covariate Xj with the treatment D (after partialling out the effect of the remaining covariates X, excluding Xj).
- `r2yxj_dx`: float with the partial R2 of covariate Xj with the outcome Y (after partialling out the effect of the remaining covariates X, excluding Xj).
- `benchmark_covariates`: a string or vector of strings with names of the variables to use for benchmark bounding.
- `kd` (default: 1): a float or list of floats with each being a multiple of the strength of association between a benchmark variable and the treatment variable to test with benchmark bounding.
- `ky` (default: nothing): same as kd except measured in terms of strength of association with the outcome variable.
"""
function ovb_partial_r2_bound(; model::Union{StatsModels.TableRegressionModel, Nothing} = nothing, treatment::Union{Nothing, String} = nothing, 
    r2dxj_x::Union{Vector{<:Real}, Real, Nothing} = nothing, r2yxj_dx::Union{Vector{<:Real}, Real, Nothing} = nothing, 
    benchmark_covariates::Union{Nothing, String, Array{String}, Dict{String, String}} = nothing, kd::Union{Vector{<:Real}, Real} = 1, 
    ky::Union{Vector{<:Real}, Real, Nothing} = nothing)

    if (isnothing(model) || isnothing(treatment)) && (isnothing(r2dxj_x) || isnothing(r2yxj_dx))
        throw(ArgumentError("ovb_partial_r2_bound requires either a GLM LinearModel objecto and a treatment name or the partial R2 values with the benchmark covariate, r2dxj_x and r2yz_dx"))
    end
    if isnothing(benchmark_covariates) && !isnothing(r2dxj_x)
        benchmark_covariates = ["manual"]
    end
    if isnothing(benchmark_covariates)
        return nothing
    elseif benchmark_covariates isa String
        benchmark_covariates = [benchmark_covariates]
    end

    if !isnothing(model)
        m = DataFrame(model.mm.m, Symbol.(coefnames(model)))
        formula_d_x = Term(Symbol(treatment)) ~ sum(Term.(Symbol.(names(m[:, Not(treatment)]))))
        treatment_fit = lm(formula_d_x, m)

        r2dxj_x = []
        r2yxj_dx = []
        if benchmark_covariates isa Array
            for b in benchmark_covariates
                push!(r2dxj_x, group_partial_r2(model = treatment_fit, covariates = [b])[1])
                push!(r2yxj_dx, group_partial_r2(model = model, covariates = [b])[1])
            end
        elseif benchmark_covariates isa Dict
            for b in keys(benchmark_covariates)
                push!(r2dxj_x, group_partial_r2(model = treatment_fit, covariates = [benchmark_covariates[b]])[1])
                push!(r2yxj_dx, group_partial_r2(model = model, covariates = [benchmark_covariates[b]])[1])
            end
        end
    elseif !isnothing(r2dxj_x)
        if r2dxj_x isa Real
            r2dxj_x = [r2dxj_x]
        end
        if r2yxj_dx isa Real
            r2yxj_dx = [r2yxj_dx]
        end
    end

    bounds = DataFrame()
    for i in eachindex(benchmark_covariates)
        r2dxj_x[i], r2yxj_dx[i] = check_r2(r2dxj_x[i], r2yxj_dx[i])
        if isnothing(ky)
            ky = kd
        end
        r2dz_x = kd * (r2dxj_x[i] / (1 - r2dxj_x[i]))
        if count(i -> (i >= 1), r2dz_x) > 0
            throw(ErrorException("Implied bound on r2dz_x >= 1, Impossible kd value. Try a lower kd"))
        end
        r2zxj_xd = kd * (r2dxj_x[i] ^ 2) ./ ((1 .- kd * r2dxj_x[i]) * (1 - r2dxj_x[i]))
        if count(i -> (i >= 1), r2zxj_xd) > 0
            throw(ErrorException("Impossible kd value. True a lower kd"))
        end
        r2yz_dx = ((sqrt.(ky) .+ sqrt.(r2zxj_xd)) ./ sqrt.(1 .- r2zxj_xd)) .^ 2 * (r2yxj_dx[i] / (1 - r2yxj_dx[i]))
        if count(i -> (i >= 1), r2yz_dx) > 0
            print("Warning: Implied bound on r2yz_dx greater than 1, try lower kd and/or ky. Setting r2yz.dx to 1.")
            if r2yz_dx isa Array
                r2yz_dx[r2yz_dx .>= 1] .= 1
            else
                if r2yz_dx >= 1
                    r2yz_dx = 1
                end
            end
        end
        if kd isa Real
            bound_label = label_maker(benchmark_covariate = benchmark_covariates[i], kd = kd, ky = ky)
            append!(bounds, Dict("bound_label" => bound_label, "r2dz_x" => r2dz_x, "r2yz_dx" => r2yz_dx))
        else
            for j in eachindex(kd)
                bound_label = label_maker(benchmark_covariate = benchmark_covariates[i], kd = kd[j], ky = ky[j])
                append!(bounds, Dict("bound_label" => bound_label, "r2dz_x" => r2dz_x[j], "r2yz_dx" => r2yz_dx[j]))
            end
        end
    end

    return bounds
end

function label_maker(; benchmark_covariate, kd, ky, digits = 2)

    if isnothing(benchmark_covariate)
        return "manual"
    else
        variable_text = join([" ", benchmark_covariate])
    end
    if ky == kd
        multiplier_text = string(round(ky, digits = digits))
    else
        multiplier_text = join([string(round(kd, digits = digits)), "/", string(round(ky, digits = digits))])
    end
    bound_label = join([multiplier_text, "x", variable_text])
    return bound_label
end

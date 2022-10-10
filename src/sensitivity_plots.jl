function ovb_contour_plot(estimate::Float64, se::Float64, dof::Int64; r2dz_x::Union{Array{<:Real}, Real, Nothing} = nothing, 
    r2yz_dx::Union{Array{<:Real}, Real, Nothing} = nothing, sensitivity_of::String = "estimate", kd::Union{Array, Real} = 1, 
    ky::Union{Array, Real, Nothing} = nothing, benchmark_covariates = nothing, bound_label = nothing, reduce::Bool = true, estimate_threshold = 0, t_threshold = 2, 
    lim = nothing, lim_y = nothing, col_contour = "black", col_thr_line = "red", label_text::Bool = true, label_bump_x = nothing, label_bump_y = nothing, 
    xlab = nothing, ylab = nothing, plot_margin_fraction = 0.05, round_dig = 3, n_levels = nothing)

    if !(sensitivity_of in ["estimate", "t-value"])
        throw(ArgumentError("sensitivity_of option must be either \"estimate\" or \"t-value\""))
    end
    estimate, r2dz_x, r2yz_dx, lim, lim_y, label_bump_x, label_bump_y = check_params(estimate, r2dz_x, r2yz_dx, lim, lim_y, label_bump_x, label_bump_y)
    plot_env = Dict()
    plot_env["lim"] = lim
    plot_env["lim_y"] = lim_y
    plot_env["reduce"] = reduce
    plot_env["sensitivity_of"] = sensitivity_of
    # add plot_env["treatment"] in another method
    
    grid_values_x = collect(0:lim/400:lim)
    grid_values_y = collect(0:lim_y/400:lim_y)
    bound_value = nothing

    z_axis = []
    if sensitivity_of == "estimate"
        for i in grid_values_x
            for j in grid_values_y
                push!(z_axis, adjusted_estimate(i, j, estimate = estimate, se = se, dof = dof, reduce = reduce)[1])
            end
        end
        threshold = estimate_threshold
        plot_estimate = estimate
        if !isnothing(r2dz_x)
            bound_value = adjusted_estimate(r2dz_x, r2yz_dx, estimate = estimate, se = se, dof = dof, reduce = reduce)
        end
    else
        for i in grid_values_x
            for j in grid_values_y
                push!(z_axis, adjusted_t(i, j, estimate = estimate, se = se, dof = dof, reduce = reduce, h0 = estimate_threshold)[1])
            end
        end
        threshold = t_threshold
        plot_estimate = (estimate - estimate_threshold) / se
        if !isnothing(r2dz_x)
            bound_value = adjusted_t(r2dz_x, r2yz_dx, estimate = estimate, se = se, dof = dof, reduce = reduce, h0 = estimate_threshold)
        end
    end
    if (bound_value isa Array) && (length(bound_value) == 1)
        bound_value = bound_value[1]
    end

    # z_axis = reshape(z_axis, length(grid_values_x), length(grid_values_y))
    z_axis = reshape(z_axis, length(grid_values_y), length(grid_values_x))
    fig, ax = subplots(1, 1, figsize = (6, 6))
    if !isnothing(n_levels)
        n_levels = n_levels - 1
    end
    CS = ax.contour(grid_values_x, grid_values_y, z_axis, colors = col_contour, linewidths = 1.0, linestyles = "solid", levels = n_levels)
    round_thr = round(threshold, digits = 2)
    if any(isapprox.(round_thr, CS.levels, atol = 0.1))
        threshold_index = findall(isapprox.(round_thr, CS.levels, atol = 0.1))[1]
        popat!(CS.collections, threshold_index)
        ax.clabel(CS, inline = 1, fontsize = 8, fmt = "%1.3g", colors = "gray", levels = CS.levels[Not(threshold_index)])
    else
        ax.clabel(CS, inline = 1, fontsize = 8, fmt = "%1.3g", colors = "gray", levels = CS.levels)
    end

    CS = ax.contour(grid_values_x, grid_values_y, z_axis, colors = col_thr_line, linewidths = 1.0, linestyles = [(0, (7, 3))], levels = [threshold])
    ax.clabel(CS, inline = 1, fontsize = 8, fmt = "%1.3g", colors = "gray")

    ax.scatter([0], [0], c = "k", marker = "^")
    ax.annotate(join(["Unadjusted\n", string(round(plot_estimate, digits = 3))]), (0.0 + label_bump_x, 0.0 + label_bump_y))

    if isnothing(xlab)
        xlab = L"Partial $R^2$ of confounder(s) with the treatment"
    end
    if isnothing(ylab)
        ylab = L"Partial $R^2$ of confounder(s) with the outcome"
    end
    xlabel(xlab)
    ylabel(ylab)
    xlim(-(lim / 15), lim)
    ylim(-(lim_y / 15), lim_y)

    if !isnothing(r2dz_x)
        r2dz_x, r2yz_dx = check_r2(r2dz_x, r2yz_dx)
        if kd isa Real
            kd = [kd]
        end
        if isnothing(ky)
            ky = kd
        end
        if isnothing(bound_label)
            bound_label = []
            for i in eachindex(kd)
                push!(bound_label, label_maker(benchmark_covariate = benchmark_covariates, kd = kd[i], ky = ky[i]))
            end
        end
        if r2dz_x isa Real
            push!(bound_label, label_maker(benchmark_covariate = nothing, kd = 1, ky = 1))
        elseif length(r2dz_x) > length(kd)
            len_dif = length(r2dz_x) - length(kd)
            for i in 1:len_dif
                push!(bound_label, label_maker(benchmark_covariate = nothing, kd = 1, ky = 1))
            end
        end
        add_bound_to_contour(r2dz_x = r2dz_x, r2yz_dx = r2yz_dx, bound_value = bound_value, bound_label = bound_label, sensitivity_of = sensitivity_of, 
        label_text = label_text, label_bump_x = label_bump_x, label_bump_y = label_bump_y, round_dig = round_dig)
    end

    x_plot_margin = plot_margin_fraction * lim
    y_plot_margin = plot_margin_fraction * lim_y
    x0, x1, y0, y1 = axis()
    axis((x0, x1 + x_plot_margin, y0, y1 + y_plot_margin))
    tight_layout()

end

function ovb_contour_plot(model::StatsModels.TableRegressionModel, treatment::String; r2dz_x = nothing, r2yz_dx = nothing, sensitivity_of::String = "estimate", 
    kd::Union{Vector, Real} = 1, ky::Union{Vector, Real, Nothing} = nothing, benchmark_covariates = nothing, bound_label = nothing, reduce::Bool = true, 
    estimate_threshold = 0, t_threshold = 2, lim = nothing, lim_y = nothing, col_contour = "black", col_thr_line = "red", label_text::Bool = true, 
    label_bump_x = nothing, label_bump_y = nothing, xlab = nothing, ylab = nothing, plot_margin_fraction = 0.05, round_dig = 3, n_levels = nothing)
    
    estimate, se, dof, r2dz_x, r2yz_dx = extract_from_model(model, treatment, benchmark_covariates, kd, ky, r2dz_x, r2yz_dx)

    ovb_contour_plot(estimate, se, dof, r2dz_x = r2dz_x, r2yz_dx = r2yz_dx, sensitivity_of = sensitivity_of, kd = kd, ky = ky, 
    benchmark_covariates = benchmark_covariates, bound_label = bound_label, reduce = reduce, estimate_threshold = estimate_threshold, t_threshold = t_threshold, 
    lim = lim, lim_y = lim_y, col_contour = col_contour, col_thr_line = col_thr_line, label_text = label_text, label_bump_x = label_bump_x, label_bump_y = label_bump_y, 
    xlab = xlab, ylab = ylab, plot_margin_fraction = plot_margin_fraction, round_dig = round_dig, n_levels = n_levels)

end

function extract_from_model(model, treatment, benchmark_covariates, kd, ky, r2dz_x, r2yz_dx)

    if isnothing(ky)
        ky = kd
    end
    if length(kd) != length(ky)
        throw(ArgumentError("kd and ky must be the same length"))
    end

    model_data = model_helper(model, treatment)
    estimate = model_data["estimate"][1]
    se = model_data["se"][1]
    dof = model_data["dof"][1]

    if isnothing(benchmark_covariates)
        return estimate, se, dof, r2dz_x, r2yz_dx
    else
        bench_bounds = obv_bounds(model, treatment, benchmark_covariates = benchmark_covariates, kd = kd, ky = ky, adjusted_estimates = false)
        if isnothing(r2dz_x)
            bounds = bench_bounds
        else
            if isnothing(r2yz_dx)
                r2yz_dx = r2dz_x
            end
            if r2dz_x isa Real
                bounds = DataFrame(Dict("r2dz_x" => [r2dz_x], "r2yz_dx" => [r2yz_dx]))
                bounds = vcat(bench_bounds, bounds, cols = :union)
            else
                bounds = DataFrame(Dict("r2dz_x" => r2dz_x, "r2yz_dx" => r2yz_dx))
                bounds = vcat(bench_bounds, bounds, cols = :union)
            end
        end
    end
    return estimate, se, dof, bounds[:, "r2dz_x"], bounds[:, "r2yz_dx"]
    
end

function extract_from_sense_obj(sense_obj::Senseobj)

    treatment = sense_obj.treatment
    estimate = sense_obj.estimate
    q = sense_obj.q
    reduce = sense_obj.reduce
    alpha = sense_obj.alpha
    se = sense_obj.se
    dof = sense_obj.dof
    benchmark_covariates = sense_obj.benchmark_covariates
    kd = sense_obj.kd
    ky = sense_obj.ky
    if reduce
        thr = estimate * (1 - q)
    else
        thr = estimate * (1 + q)
    end
    t_thr = abs(quantile(TDist(dof - 1), alpha / 2)) * sign(sense_obj.sensitivity_statistics["t_statistic"][1])

    if isnothing(sense_obj.bounds)
        r2dz_x = nothing
        r2yz_dx = nothing
        bound_label = ""
    else
        r2dz_x = sense_obj.bounds[:, "r2dz_x"]
        r2yz_dx = sense_obj.bounds[:, "r2yz_dx"]
        bound_label = sense_obj.bounds[:, "bound_label"]
    end

    return treatment, estimate, se, dof, r2dz_x, r2yz_dx, bound_label, reduce, thr, t_thr, benchmark_covariates, kd, ky
end

function add_bound_to_contour(; kd = 1, ky = nothing, r2dz_x, r2yz_dx, bound_value = nothing, bound_label = nothing, sensitivity_of = nothing, 
    label_text = true, label_bump_x = nothing, label_bump_y = nothing, round_dig = 3)

    if r2dz_x isa Real
        r2dz_x = [r2dz_x]
    end
    if r2yz_dx isa Real
        r2yz_dx = [r2yz_dx]
    end
    if bound_value isa Real
        bound_value = [bound_value]
    end

    for i in eachindex(r2dz_x)
        scatter(r2dz_x[i], r2yz_dx[i], c = "red", marker = "D", edgecolors = "black")
        if label_text
            if bound_label isa Real
                bound_label = [bound_label]
            end
            if !isnothing(bound_label) && !isnothing(bound_value)
                bound_value[i] = round(bound_value[i], digits = round_dig)
                label = join([string(bound_label[i]), "\n(", string(bound_value[i]), ")"])
            else
                label = bound_label[i]
            end
            annotate(label, (r2dz_x[i] + label_bump_x, r2yz_dx[i] + label_bump_y))
        end
    end

end

function check_params(estimate, r2dz_x, r2yz_dx, lim, lim_y, label_bump_x, label_bump_y)

    if isnothing(r2yz_dx)
        r2yz_dx = r2dz_x
    end

    if isnothing(lim)
        if isnothing(r2dz_x)
            lim = 0.4
        else
            if r2yz_dx isa Real
                lim = minimum([maximum([r2dz_x * 1.2, 0.4]), 1 - 1e-12])
            else
                lim = minimum([maximum(vcat(r2dz_x * 1.2, 0.4)), 1 - 1e-12])
            end
        end
    end
    if isnothing(lim_y)
        if isnothing(r2yz_dx)
            lim_y = 0.4
        else
            if r2yz_dx isa Real
                lim_y = minimum([maximum([r2yz_dx * 1.2, 0.4]), 1 - 1e-12])
            else
                lim_y = minimum([maximum(vcat(r2yz_dx * 1.2, 0.4)), 1 - 1e-12])
            end
        end
    end
    if isnothing(label_bump_x)
        label_bump_x = lim / 30
    end
    if isnothing(label_bump_y)
        label_bump_y = lim_y / 30
    end
    if lim > 1
        lim = 1 - 1e-12
        print("Warning: contour limit larger than 1 was set to 1")
    elseif lim < 0
        lim = 0.4
        print("Warning: contour limit lesser than 0 was set to 0")
    end
    if lim_y > 1
        lim_y = 1 - 1e-12
        print("Warning: contour limit larger than 1 was set to 1")
    elseif lim_y < 0
        lim_y = 0.4
        print("Warning: contour limit lesser than 0 was set to 0")
    end
    
    return estimate, r2dz_x, r2yz_dx, lim, lim_y, label_bump_x, label_bump_y
end

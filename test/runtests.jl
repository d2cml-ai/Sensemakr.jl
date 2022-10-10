# cd("/Documents/Work/d2cml-ai/Sensemakr.jl")

using Test
using Sensemakr
using CSV, DataFrames, GLM

const path = joinpath(dirname(@__FILE__), "..", "data", "darfur.csv");

@testset "Sensemakr.jl" begin
    # set up data
    darfur = CSV.read(path, DataFrame);
    form = @formula(peacefactor ~ directlyharmed + age + farmer_dar + herder_dar + pastvoted + hhsize_darfur + female + village);
    fitted_model = lm(form, darfur);
    atol = 0.000001;

    #tests
    @test robustness_value(model = fitted_model, covariates = "directlyharmed")[1] ≈ 0.138776 atol = atol
    @test robustness_value(model = fitted_model, covariates = "directlyharmed", q = 1/2)[1] ≈ 0.072027 atol = atol = atol
    @test robustness_value(model = fitted_model, covariates = "directlyharmed", q = 1/2, alpha = 0.05)[1] ≈ 0.004563 atol = atol
    @test robustness_value(t_statistic = 4.18445, dof = 783)[1] ≈ 0.138776 atol = atol

    @test partial_r2(model = fitted_model, covariates = "directlyharmed")[1] ≈ 0.021873093341109 atol = atol
    @test partial_r2(model = fitted_model, covariates = "female")[1] ≈ 0.10903391542784427 atol = atol
    @test partial_r2(t_statistic = 4.18445, dof = 783) ≈ 0.021873093496457607 atol = atol

    @test partial_f2(model = fitted_model, covariates = "directlyharmed")[1] ≈ 0.022362 atol = atol
    @test partial_f2(model = fitted_model, covariates = "female")[1] ≈ 0.122377 atol = atol
    @test partial_f2(t_statistic = 4.18445, dof = 783) ≈ 0.022362224524265645 atol = atol

    @test group_partial_r2(model = fitted_model, covariates = ["female", "pastvoted"])[1] ≈ 0.11681276064557611 atol = atol
    
    @test sensitivity_stats(model = fitted_model, treatment = "directlyharmed")["dof"] == 783

    @test adjusted_estimate(0.05, 0.05, model = fitted_model, treatment = "directlyharmed")[1] ≈ 0.06393214421078033 atol = atol
    
    @test adjusted_se(0.05, 0.05, model = fitted_model, treatment = "directlyharmed")[1] ≈ 0.02327140296812648 atol = atol

    @test adjusted_t(0.05, 0.05, model = fitted_model, treatment = "directlyharmed")[1] ≈ 2.7472406497513093 atol = atol

    @test obv_bounds(fitted_model, "directlyharmed", benchmark_covariates = ["female", "pastvoted"], kd = [1, 2, 3]) isa DataFrame
end


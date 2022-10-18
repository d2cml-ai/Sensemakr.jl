var documenterSearchIndex = {"docs":
[{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"EditURL = \"<unknown>/quickstart.jl\"","category":"page"},{"location":"quickstart/#Quickstart","page":"Quickstart","title":"Quickstart","text":"","category":"section"},{"location":"quickstart/#Introduction","page":"Quickstart","title":"Introduction","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"In this quickstart, we demonstrate the basic usage of the package by reproducing Section 5 of Cinelli and Hazlett (2020), which estimates the effects of exposure to violence on attitudes towards peace, in Darfur. Throughout this manual, we mainly focus on the code. For detailed explanations, please refer to the original paper.","category":"page"},{"location":"quickstart/#Violence-in-Darfur","page":"Quickstart","title":"Violence in Darfur","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"In 2003 and 2004, the Darfurian government orchestrated a horrific campaign of violence against civilians, killing an estimated two hundred thousand people. In this application, we are interested in answering the following question: did being directly exposed to harm make individuals more “angry,” and thus more likely to ask for revenge, or did it make them more \"weary,\" and thus more likely to ask for peace?","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"Sensemakr comes with the Darfur dataset, which can be loaded with the command data.load_darfur(). More details about the data can be found in the documentation, as well as in Hazlett (2019) and Cinelli and Hazlett (2020).","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"# imports Sensemakr\nusing Sensemakr\n# loads darfur data\nusing DataFrames\ndarfur = load_darfur();\nfirst(darfur, 5)","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"<div class=\"data-frame\"><p>5 rows × 14 columns</p><table class=\"data-frame\"><thead><tr><th></th><th>wouldvote</th><th>peacefactor</th><th>peace_formerenemies</th><th>peace_jjindiv</th><th>peace_jjtribes</th><th>gos_soldier_execute</th><th>directlyharmed</th><th>age</th><th>farmer_dar</th><th>herder_dar</th><th>pastvoted</th><th>hhsize_darfur</th><th>village</th><th>female</th></tr><tr><th></th><th title=\"Int64\">Int64</th><th title=\"Float64\">Float64</th><th title=\"Int64\">Int64</th><th title=\"Int64\">Int64</th><th title=\"Int64\">Int64</th><th title=\"Int64\">Int64</th><th title=\"Int64\">Int64</th><th title=\"Int64\">Int64</th><th title=\"Int64\">Int64</th><th title=\"Int64\">Int64</th><th title=\"Int64\">Int64</th><th title=\"Int64\">Int64</th><th title=\"InlineStrings.String31\">String31</th><th title=\"Int64\">Int64</th></tr></thead><tbody><tr><th>1</th><td>0</td><td>1.0</td><td>1</td><td>1</td><td>1</td><td>0</td><td>0</td><td>30</td><td>0</td><td>0</td><td>1</td><td>23</td><td>Abdel Khair</td><td>0</td></tr><tr><th>2</th><td>0</td><td>0.706831</td><td>0</td><td>1</td><td>1</td><td>0</td><td>0</td><td>20</td><td>1</td><td>0</td><td>1</td><td>5</td><td>Abdi Dar</td><td>1</td></tr><tr><th>3</th><td>1</td><td>0.0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>0</td><td>45</td><td>1</td><td>0</td><td>0</td><td>15</td><td>Abu Sorog</td><td>0</td></tr><tr><th>4</th><td>1</td><td>0.495178</td><td>1</td><td>0</td><td>0</td><td>0</td><td>1</td><td>55</td><td>0</td><td>0</td><td>0</td><td>9</td><td>Abu Dejaj</td><td>0</td></tr><tr><th>5</th><td>0</td><td>0.0</td><td>0</td><td>0</td><td>0</td><td>1</td><td>1</td><td>25</td><td>1</td><td>0</td><td>1</td><td>7</td><td>Abu Dejaj</td><td>1</td></tr></tbody></table></div>","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"A natural approach for such problem is to run the following linear regression model, where we regress peacefactor on directlyharmed, further adjusting for village, female as well as other covariates. Here we run this regression using GLM","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"using GLM\nform = @formula(peacefactor ~ directlyharmed + age + farmer_dar + herder_dar + pastvoted + hhsize_darfur + female + village);\nfitted_model = lm(form, darfur);","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"The above regression results in the following estimate and standard errors for the coefficient of directlyharmed:","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"first(coeftable(fitted_model), 2)[2]","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"(Name = \"directlyharmed\", var\"Coef.\" = 0.09731581928495726, var\"Std. Error\" = 0.023256537809811153, t = 4.184449984808271, var\"Pr(>|t|)\" = 3.182339959789431e-5, var\"Lower 95%\" = 0.05166327477010026, var\"Upper 95%\" = 0.14296836379981426)","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"According to this model, those individual who were directly exposed to harm became on average more “pro-peace,” not less.","category":"page"},{"location":"quickstart/#Sensitivity-Analysis","page":"Quickstart","title":"Sensitivity Analysis","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"The causal interpretation of the previous estimate, however, relies on the assumption that village and gender are sufficient for control of confounding—in other words, it requires the assumption of no unobserved confounders. What if is wrong? How strong would these unobserved variables have to be in order to change the original research conclusions? The goal of Sensemakr is precisely that, i.e, to make it easier to understand the impact that omitted variables would have on a regression result.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"The main function of the package is the function sensemakr. The sensemakr function defines an object of type sensemakr, performing the most commonly required sensitivity analyses which can then be further explored with the summary and plot methods of the object. This function is mainly a convenience wrapper for other sensitivity functions defined in the package, which can also be called directly, as we detail later in the documentation.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"In the code chunk below, we apply the function sensemakr to our OLS model from statsmodel.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"darfur_sense = sensemakr(fitted_model, \"directlyharmed\", benchmark_covariates = \"female\", kd = [1, 2, 3], ky = [1, 2, 3], q = 1.0, alpha = 0.05, reduce = true);","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"The main arguments of the call are:","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"model: the OLS model with the outcome regression. In our case, darfur_model.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"treatment: the name of the treatment variable. In our case, \"directlyharmed\".","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"benchmark_covariates: the names of covariates that will be used to bound the plausible strength of the unobserved confounders. Here, we put \"female\", which is arguably one of the main determinants of exposure to violence, and also a strong determinant of attitudes towards peace.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"kd and ky: these arguments parameterize how many times stronger the confounder is related to the treatment (kd) and to the outcome (ky) in comparison to the observed benchmark covariates (in this case, female). In our example, setting kd = [1, 2, 3] and ky = [1, 2, 3] means we want to investigate the maximum strength of a confounder once, twice, or three times as strong as female (in explaining treatment and outcome variation). If only kd is given, ky will be set equal to kd.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"q: fraction of the effect estimate that would have to be explained away to be problematic. Setting q = 1, means that a reduction of 100% of the current effect estimate, that is, a true effect of zero, would be deemed problematic. The default is 1.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"alpha: significance level of interest for statistical inference. The default is 0.05.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"reduce: should we consider confounders acting towards increasing or reducing the absolute value of the estimate? The default is reduce = true, which means we are considering confounders that pull the estimate towards (or through) zero.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"Using the default arguments, one can simplify the previous call to:","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"darfur_sense = sensemakr(fitted_model, \"directlyharmed\", benchmark_covariates = \"female\", kd = [1, 2, 3]);","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"Once we run sensemakr, we can now explore the sensitivity analysis results.","category":"page"},{"location":"quickstart/#Minimal-sensitivity-reporting","page":"Quickstart","title":"Minimal sensitivity reporting","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"The print method of Sensemakr provides a quick review of the original (unadjusted) estimates along with three summary sensitivity statistics suited for routine reporting: the partial R^2 of the treatment with the outcome, the robustness value (RV) required to reduce the estimate entirely to zero (i.e. q=1), and the RV beyond which the estimate would no longer be statistically distinguishable from zero at the 0.05 level","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"print(darfur_sense)","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"Sensitivity Analysis to Unobserved Confounding\n\nModel Formula: hhsize_darfur ~ farmer_dar + herder_dar + female + village + peacefactor + age + directlyharmed + pastvoted\n\nNull hypothesis: q = 1.0 and reduce = true\n\nUnadjusted Estimates of \"directlyharmed\":\n   Coef. Estimate: 0.097\n   Standard Error: 0.023\n   t-value: 4.184\n\nSensitivity Statistics:\n   Partial R2 of treatment with outcome: 0.022\n   Robustness Value, q = 1.0: 0.139\n   Robustness Value, q = 1.0 alpha = 0.05: 0.076\n\n","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"The package also provides a function that outputs code for a latex or html table with these results. If used in an interactive environment, such as a Jupyter notebook, the table is also automatically displayed in the notebook.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"# html code for minimal reporting table","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"These three sensitivity statistics provide a minimal reporting for sensitivity analysis. More precisely:","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"The robustness value for bringing the point estimate of directlyharmed exactly to zero (RV_q=1) is 13.9% . This means that unobserved confounders that explain 13.9% of the residual variance both of the treatment and of the outcome are sufficiently strong to explain away all the observed effect. On the other hand, unobserved confounders that do not explain at least 13.9% of the residual variance both of the treatment and of the outcome are not sufficiently strong to do so.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"The robustness value for testing the null hypothesis that the coefficient of directlyharmed is zero (RV_q=1 alpha = 005) falls to 7.6%. This means that unobserved confounders that explain 7.6% of the residual variance both of the treatment and of the outcome are sufficiently strong to bring the lower bound of the confidence interval to zero (at the chosen significance level of 5%). On the other hand, unobserved confounders that do not explain at least 7.6% of the residual variance both of the treatment and of the outcome are not sufficiently strong to do so.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"Finally, the partial R^2 of directlyharmed with peacefactor means that, in an extreme scenario, in which we assume that unobserved confounders explain all of the left out variance of the outcome, these unobserved confounders would need to explain at least 2.2% of the residual variance of the treatment to fully explain away the observed effect.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"The lower corner of the table, further provides bounds on the strength of an unobserved confounder as strong as the observed covariate female, resulting in R^2_Ysim Zmid X D = 125 and R^2_Dsim Zmid X = 9. Since both of those are below the RV of 13.9%, we conclude confounders as strong as female, in explaining treatment and outcome variations, are not sufficiently strong to explain away the observed estimate.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"Moreover, the bound of R^2_Dsim Zmid X = 9 is below the partial 𝑅2 of the treatment with the outcome, R^2_Ysim Dmid X = 22, this means that even an extreme confounder explaining all residual variation of the outcome, and as strongly associated with the treatment as female would not be able to overturn the research conclusions.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"The summary function of Sensemakr produces verbose output similar to the text explanations above, so that researchers can directly cite or include such texts in their reports.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"summary(darfur_sense)","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"Sensitivity Analysis to Unobserved Confounding\n\nModel Formula: hhsize_darfur ~ farmer_dar + herder_dar + female + village + peacefactor + age + directlyharmed + pastvoted\n\nNull hypothesis: q = 1.0 and reduce = true\n-- This means we are considering biases that reduce the absolute value of the current estimate\n-- The null hypothesis deemed problematic is H0:tau = 0.0\n\nUnadjusted Estimates of \"directlyharmed\":\n   Coef. Estimate: 0.097\n   Standard Error: 0.023\n   t-value: 4.184\nSensitivity Statistics:\n   Partial R2 of treatment with outcome: 0.022\n   Robustness Value, q = 1.0: 0.139\n   Robustness Value, q = 1.0 alpha = 0.05: 0.076\n\nVerbal interpretation of sensitivity statistics:\n\n-- Partial R2 of the treatment with the outcome: an extreme confounder (orthogonal to the covariates) that explains 100% of the residual variance of the outcome, would need to explain at least 2.187 % of the residual variance of the treatment to fully account for the observed estimated effect.\n\n-- Robustness Value, q = 1.0: unobserved confounders (orthogonal to the covariates) that of both the treatment and the outcome are strong enough to bring the point estimate to 0.0 (a bias of 100.0% of the original estimate). Conversely, unobserved confounders that do not explain more than 13.878% of the residual variance of both the treatment and the outcome are not strong enough to bring the point estimate to 0.0.\n\n-- Robustness Value,q = 1.0, alpha = 0.05: unobserved confounders (orthogonal to the covariates) that explain more than 7.626% of the residual variance of both the treatment and the outcome are strong enough to bring the estimate to a range where it is no longer 'statistically different' from 0.0 (a bias of 100.0% of the original estimate), at the significance level of alpha = 0.05. Conversely, unobserved confounders that do not explain more than7.626% of the residual varianceof both the treatment and the outcome are not strong enough to bring the estimate to a range where it is no longer 'statistically different' from 0.0, at the significance level of alpha = 0.05.\n\nBounds on omitted variable bias:\n--The table below shows the maximum strength of unobserved confounders with association with the treatment and the outcome bounded by a multiple of the observed explanatory power of the chosen benchmark covariate(s).\n\n3×9 DataFrame\n Row │ bound_label  r2dz_x      r2yz_dx   treatment       adjusted_estimate  adjusted_se  adjusted_t  adjusted_lower_CI  adjusted_upper_CI\n     │ String       Float64     Float64   String          Float64            Float64      Float64     Float64            Float64\n─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n   1 │ 1.0x female  0.00916429  0.124641  directlyharmed          0.0752203    0.0218733     3.4389          0.032283            0.118158\n   2 │ 2.0x female  0.0183286   0.249324  directlyharmed          0.0529152    0.0203501     2.60025         0.012968            0.0928623\n   3 │ 3.0x female  0.0274929   0.37405   directlyharmed          0.030396     0.0186701     1.62806        -0.00625328          0.0670453\n","category":"page"},{"location":"quickstart/#Sensitivity-plot","page":"Quickstart","title":"Sensitivity plot","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"Using the plot method for sensemakr, we can further refine our sensitivity analysis by visually exploring the whole range of possible estimates that confounders with different strengths could cause.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"Let us begin by examining contour plots for the point estimate.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"plot(darfur_sense)","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"(Image: Figure_1)","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"The horizontal axis shows the hypothetical residual share of variation of the treatment that unobserved confounding explains, R^2_Dsim Z bf X . The vertical axis shows the hypothetical partial R^2 of unobserved confounding with the outcome, R^2_Ysim Z bf X D. The contours show what would be the estimate for directlyharmed that one would have obtained in the full regression model including unobserved confounders with such hypothetical strengths. Note the plot is parameterized in way that hurts our preferred hypothesis, by pulling the estimate towards zero—the direction of the bias was set in the argument reduce = true of sensemakr().","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"The bounds on the strength of confounding, determined by the parameter kd = [1, 2, 3] in the call for sensemakr(), are also shown in the plot. Note that the plot reveals that the direction of the effect (positive) is robust to confounding once, twice or even three times as strong as the observed covariate female, although in this last case the magnitude of the effect is reduced to a third of the original estimate.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"We now examine the sensitivity of the t-value for testing the null hypothesis of zero effect. For this, it suffices to change the option sensitivity_of = \"t-value\".","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"plot(darfur_sense, sensitivity_of = \"t-value\")","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"(Image: Figure_2)","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"The plot reveals that, at the 5% significance level, the null hypothesis of zero effect would still be rejected given confounders once or twice as strong as female. However, by contrast to the point-estimate, accounting for sampling uncertainty now means that the null hypothesis of zero effect would not be rejected with the inclusion of a confounder three times as strong as female.","category":"page"},{"location":"quickstart/#Sensitivity-to-extreme-scenarios","page":"Quickstart","title":"Sensitivity to extreme scenarios","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"Sometimes researchers may be better equipped to make plausibility judgments about the strength of determinants of the treatment assignment mechanism, and have less knowledge about the determinants of the outcome. In those cases, sensitivity plots using extreme scenarios are a useful option. These are produced with the option plot_type = \"extreme\". Here one assumes confounding explains all or some large fraction of the residual variance of the outcome, then vary how strongly such confounding is hypothetically related to the treatment, to see how this affects the resulting point estimate.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"plot(darfur_sense, plot_type = \"extreme\")","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"(Image: Figure_3)","category":"page"},{"location":"quickstart/#References","page":"Quickstart","title":"References","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"Cinelli, C. Hazlett, C. (2020) “Making Sense of Sensitivity: Extending Omitted Variable Bias”. Journal of the Royal Statistical Society, Series B (Statistical Methodology).","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"Hazlett, C. (2019). Angry or Weary? How Violence Impacts Attitudes toward Peace among Darfurian Refugees. Journal of Conflict Resolution.","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"This page was generated using Literate.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Sensemakr","category":"page"},{"location":"#Sensemakr","page":"Home","title":"Sensemakr","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Sensemakr, which implements a suite of sensitivity analysis tools that makes it easier to understand the impact of omitted variables in linear regression models, as discussed in Cinelli and Hazlett (2020).","category":"page"},{"location":"","page":"Home","title":"Home","text":"The Python version of the package closely mirrors the R version, which can be found here.","category":"page"},{"location":"#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"quickstart.md\"\n]\nDepth = 2","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Sensemakr]","category":"page"}]
}

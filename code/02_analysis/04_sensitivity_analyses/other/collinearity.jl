
cd("../../../../")

using Pkg

Pkg.activate("./")

using CSV, DataFrames, MLJ, MLJDecisionTreeInterface,
    DecisionTree, Plots, StatsBase, Random, MLJBase,
    JLD2, Shapley, Shapley.ComputationalResources,
    NaNStatistics, Loess, RData, Serialization

## load in data
data = CSV.read("data/complete_data.csv", DataFrame);
data = data[:,1:43]
data = data[completecases(data),:];

folds_formatted = load_object("data/random_forests/folds_formatted.jld2");

data.holdout .= 1
data[vcat(folds_formatted[1][1], folds_formatted[1][2]), :holdout] .= 0

data.hs_class .= "not_high";
data[data.cbi .> 2.25, :hs_class] .= "high";

## target
target = Vector(data[:,:hs_class]);
target = coerce(target, OrderedFactor);
levels!(target, ["not_high", "high"]);


##---------------------------------------------------------------
## Neighborhood Scale
##---------------------------------------------------------------

##---------------------------------------------------------------
## remove metrics with high collinearity
##---------------------------------------------------------------
function safe_parse_float(x)
    if x in ("NA", "missing", "")
        return missing
    else
        return tryparse(Float64, x)
    end
end

correl = CSV.read("data/correlation_matrices/correlation_matrix_neighborhood.csv", DataFrame)
correl.value .= safe_parse_float.(correl.value)
correl[.!ismissing.(correl.value) .&& abs.(correl.value) .> 0.5,:]

## try with just mean_dens, em, and tpi out of these correlated values

## scale and center parameters
data.clust_30_scaled = (data.clust_30 .- mean(data.clust_30)) ./ std(data.clust_30);
data.mean_dens_30_scaled = (data.mean_dens_30 .- mean(data.mean_dens_30)) ./ std(data.mean_dens_30);
data.mean_area_30_scaled = (data.mean_area_30 .- mean(data.mean_area_30)) ./ std(data.mean_area_30);
data.mean_ht_30_scaled = (data.mean_ht_30 .- mean(data.mean_ht_30)) ./ std(data.mean_ht_30);
data.em_scaled = (data.em .- mean(data.em)) ./ std(data.em);
data.cwd_scaled = (data.cwd .- mean(data.cwd)) ./ std(data.cwd);
data.slope_scaled = (data.slope .- mean(data.slope)) ./ std(data.slope);
data.tpi_scaled = (data.tpi .- mean(data.tpi)) ./ std(data.tpi);
data.heat_load_scaled = (data.heat_load .- mean(data.heat_load)) ./ std(data.heat_load);
data.prev_sev_scaled = (data.prev_sev .- mean(filter(!isnan, data.prev_sev))) ./ std(filter(!isnan, data.prev_sev));
data.avg_fuel_moisture_scaled = (data.avg_fuel_moisture .- mean(data.avg_fuel_moisture)) ./ std(data.avg_fuel_moisture);
data.hdw_scaled = (data.hdw .- mean(data.hdw)) ./ std(data.hdw);

cols = [:mean_dens_30_scaled, :em_scaled,
        :cwd_scaled, :slope_scaled, :tpi_scaled, :prev_sev_scaled, :avg_fuel_moisture_scaled,
        :hdw_scaled];
features = DataFrame(fire_name = coerce(data[:,:fire_name], OrderedFactor));

## coerce to correct scitype
for i in 1:length(cols)
    features[:,cols[i]]  = coerce(float.(identity.(data[:,cols[i]])), Continuous)
end

holdout_features = copy(features[data.holdout .== 1, :])
holdout_target = copy(target[data.holdout .== 1, :])

## load classifier
RandomForest = @MLJ.load RandomForestClassifier pkg = DecisionTree

##---------------------------------------------------------------
## Fit model with tuned parameters
##---------------------------------------------------------------

forest_model = MLJDecisionTreeInterface.RandomForestClassifier(n_trees = 40,
                                                               max_depth = 30,
                                                               n_subfeatures = 2)
fm = machine(forest_model,
             features,
             target)

## fit model, takes awhile
MLJ.fit!(fm, rows = vcat(folds_formatted[1][1], folds_formatted[1][2]))

serialize("data/forest_model_colin_30.ser", fm)

##---------------------------------------------------------------
## check out of sample prediction performance and compute feature importances
##---------------------------------------------------------------
fm = deserialize("data/forest_model_colin_30.ser")

ho_list = findall(data.holdout .== 1);

yhat = MLJ.predict(fm, rows = ho_list);

auc(yhat, data[ho_list, :hs_class])
confusion_matrix(mode.(yhat), data[ho_list, :hs_class])

##---------------------------------------------------------------
## Compute shapley metrics
##---------------------------------------------------------------
using Shapley: MonteCarlo

## create subset on which to compute shapley values
Random.seed!(1)
exp_list = sample(ho_list, 5000)
explain = copy(features[exp_list, :])

## compute shapley values
shap = shapley(η -> MLJ.predict(fm, η), MonteCarlo(CPUThreads(), 64), explain)
shap2 = shap
serialize("data/random_forests/shapley/shapley_colin_30.ser", shap)

shap = deserialize("data/random_forests/shapley/shapley_colin_30.ser")

## summarize shapley values
shapley_summary = DataFrame(feature = names(features),
                            mean_not_high = Vector{Float64}(undef, ncol(features)),
                            sd_not_high = Vector{Float64}(undef, ncol(features)),
                            mean_high = Vector{Float64}(undef, ncol(features)),
                            sd_high = Vector{Float64}(undef, ncol(features)))

for i in 1:length(shap)
    shp_mat = pdf(shap[i], ["not_high", "high"])
    shapley_summary[i, [2,4]] .= vec(mean(abs.(shp_mat), dims = 1))
    shapley_summary[i, [3,5]] .= vec(std(abs.(shp_mat), dims = 1))
end

sort(shapley_summary, order(:mean_high, rev = true))

data_plot = sort(shapley_summary, order(:mean_high, rev = true))
bar(data_plot.mean_high, orientation = :h, yticks = (1:16, data_plot.feature),
    xlims = [0, 0.15], frame = :box, legend = :none, color = :black, title = "Feature importance (high)",
    yflip = true)


savefig("plots/shapley_colin_30.pdf")

CSV.write("data/random_forests/shapley/shapley_colin_30.csv", shapley_summary)





##---------------------------------------------------------------
## FOR LARGE SCALE METRICS
##
##
##
##---------------------------------------------------------------

##---------------------------------------------------------------
## Fit and tune random forest model
##---------------------------------------------------------------

correl = CSV.read("data/correlation_matrices/correlation_matrix_stand.csv", DataFrame)
correl.value .= safe_parse_float.(correl.value)
correl[.!ismissing.(correl.value) .&& abs.(correl.value) .> 0.5,:]

data.clust_180_scaled = (data.clust_180 .- mean(data.clust_180)) ./ std(data.clust_180)
data.mean_dens_180_scaled = (data.mean_dens_180 .- mean(data.mean_dens_180)) ./ std(data.mean_dens_180)
data.mean_area_180_scaled = (data.mean_area_180 .- mean(data.mean_area_180)) ./ std(data.mean_area_180)
data.mean_ht_180_scaled = (data.mean_ht_180 .- mean(data.mean_ht_180)) ./ std(data.mean_ht_180)

cols = [:mean_dens_180_scaled, :mean_area_180_scaled, :mean_ht_180_scaled, :em_scaled,
        :cwd_scaled, :slope_scaled, :tpi_scaled, :prev_sev_scaled, :avg_fuel_moisture_scaled,
        :hdw_scaled];
features = DataFrame(fire_name = coerce(data[:,:fire_name], OrderedFactor));

## coerce to correct scitype
for i in 1:length(cols)
    features[:,cols[i]]  = coerce(float.(identity.(data[:,cols[i]])), Continuous)
end

holdout_features = copy(features[data.holdout .== 1, :])
holdout_target = copy(target[data.holdout .== 1, :])

## load classifier
RandomForest = @MLJ.load RandomForestClassifier pkg = DecisionTree

##---------------------------------------------------------------
## Fit model
##---------------------------------------------------------------

forest_model_180 = MLJDecisionTreeInterface.RandomForestClassifier(n_trees = 40,
                                                               max_depth = 40,
                                                               n_subfeatures = 3)
fm_180 = machine(forest_model_180,
             features,
             target)

## fit model, takes awhile\
MLJ.fit!(fm_180, rows = vcat(folds_formatted[1][1], folds_formatted[1][2]))

serialize("data/random_forests/forest_model_colin_180.ser", fm_180)

##---------------------------------------------------------------
## check out of sample prediction performance and compute feature importances
##---------------------------------------------------------------

fm_180 = deserialize("data/random_forests/forest_model_colin_180.ser")

ho_list = findall(data.holdout .== 1)

yhat_180 = MLJ.predict(fm_180, rows = ho_list)

auc(yhat_180, data[ho_list, :hs_class])
confusion_matrix(mode.(yhat_180), data[ho_list, :hs_class])


##---------------------------------------------------------------
## Compute shapley metrics
##---------------------------------------------------------------

## create subset on which to compute shapley values
Random.seed!(2)
exp_list = sample(ho_list, 5000)
explain = copy(features[exp_list, :])

## compute shapley values
shap_180 = shapley(η -> MLJ.predict(fm_180, η), MonteCarlo(CPUThreads(), 64), explain)

serialize("data/random_forests/shapley/shapley_colin_180.ser", shap_180)

## summarize shapley values
shapley_summary_180 = DataFrame(feature = names(features),
                                mean_not_high = Vector{Float64}(undef, ncol(features)),
                                sd_not_high = Vector{Float64}(undef, ncol(features)),
                                mean_high = Vector{Float64}(undef, ncol(features)),
                                sd_high = Vector{Float64}(undef, ncol(features)))

for i in 1:length(shap_180)
    shp_mat = pdf(shap_180[i], ["not_high", "high"])
    shapley_summary_180[i, [2,4]] .= vec(mean(abs.(shp_mat), dims = 1))
    shapley_summary_180[i, [3,5]] .= vec(std(abs.(shp_mat), dims = 1))
end

sort(shapley_summary_180, order(:mean_high, rev = true))

data_plot = sort(shapley_summary_180, order(:mean_high, rev = true))
bar(data_plot.mean_high, orientation = :h, yticks = (1:16, data_plot.feature),
    xlims = [0, 0.15], frame = :box, legend = :none, color = :black, title = "Feature importance (high)",
    yflip = true)
savefig("plots/shapley_colin_180.pdf")
CSV.write("data/random_forests/shapley/shapley_colin_180.csv", shapley_summary_180)

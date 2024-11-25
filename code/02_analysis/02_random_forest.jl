
cd("../../")

using Pkg

Pkg.activate("./")

using CSV, DataFrames, MLJ, MLJDecisionTreeInterface,
    DecisionTree, Plots, StatsBase, Random, MLJBase,
    JLD2, Shapley, Shapley.ComputationalResources,
    NaNStatistics, Loess, RData

## load in data
data = CSV.read("data/complete_data.csv", DataFrame);
data = data[completecases(data),:];

##---------------------------------------------------------------
## create grid for kfold cross validation
##---------------------------------------------------------------

bb = RData.load("data/autocor_scale_m.rds")

fire_ids = unique(data.fire_name);

## initialize dataframe
boxes_full = DataFrame(lower_left = Vector{Tuple{Float64, Float64}}(undef, 0), id = Vector{Int64}(undef, 0),
                       fire_name = Vector{String}(undef, 0))

for f in fire_ids

    subdata = data[data.fire_name .== f, :]
    xmin = minimum(subdata.x)
    xmax = maximum(subdata.x)
    ymin = minimum(subdata.y)
    ymax = maximum(subdata.y)

    xc = range(xmin, xmax, step = bb)
    yc = range(ymin, ymax, step = bb)

    boxes = DataFrame(lower_left = vec(collect(Base.product(xc, yc))))
    boxes.id .= (findall(fire_ids .== f) .* 10000) .+ [1:1:nrow(boxes);]

    notempty = []
    for i in 1:nrow(boxes)
        xmn = boxes[i, :lower_left][1]
        ymn = boxes[i, :lower_left][2]
        xmx = xmn + bb
        ymx = ymn + bb

        if nrow(subdata[subdata.x .>= xmn .&& subdata.x .< xmx .&& subdata.y .>= ymn .&& subdata.y .< ymx, :]) > 0
            push!(notempty, i)
        end
    end
    boxes = boxes[notempty, :]
    boxes.fire_name .= f
    boxes_full = vcat(boxes_full, boxes)

end

plot(boxes_full.lower_left, seriestype = :scatter)
CSV.write("data/boxes_for_kfold.csv", boxes_full)


testing_prop =  0.2
training_ids = []
testing_ids = []
for f in fire_ids
    sub = boxes_full[boxes_full.fire_name .== f, :]
    Random.seed!(1)
    tr = StatsBase.sample(sub.id, Int(round((1-testing_prop) * nrow(sub))), replace = false)
    training_ids = vcat(training_ids, tr)
    for i in 1:nrow(sub)
        if sub.id[i] ∉ tr
            push!(testing_ids, sub.id[i])
        end
    end
end

## check to make sure proportions are correct
length(testing_ids) / nrow(boxes_full)
length(training_ids) / nrow(boxes_full)

## generate new dataframes with testing and training data
testing_boxes = DataFrame(lower_left = Vector{Tuple{Float64, Float64}}(undef, 0), fire_name = Vector{String}(undef, 0))
for i in 1:length(testing_ids)
    testing_boxes = vcat(testing_boxes, boxes_full[boxes_full.id .== testing_ids[i], [:lower_left, :fire_name]])
end
testing_boxes

plot!(testing_boxes.lower_left, group = testing_boxes.fire_name, seriestype = :scatter)

training_boxes = DataFrame(lower_left = Vector{Tuple{Float64, Float64}}(undef, 0),
                           fire_name = Vector{String}(undef, 0), id = Vector{Int64}(undef, 0))
for i in 1:length(training_ids)
    training_boxes = vcat(training_boxes, boxes_full[boxes_full.id .== training_ids[i], [:lower_left, :fire_name, :id]])
end
training_boxes

plot(training_boxes.lower_left, group = training_boxes.fire_name, seriestype = :scatter)

CSV.write("data/testing_boxes.csv", testing_boxes)
CSV.write("data/training_boxes.csv", training_boxes)


##---------------------------------------------------------------
## Generate folds for CV
##---------------------------------------------------------------
training_boxes = CSV.read("data/training_boxes.csv", DataFrame)
training_boxes.lower_left .= eval.(Meta.parse.(training_boxes.lower_left))

## function to generate folds for kfold CV
function gen_folds(vector, nfolds, seed)
    len = length(vector)
    n,r = divrem(len, nfolds)
    b = collect(1:n:len+1)
    for i in 1:length(b)
        b[i] += i > r ? r : i-1
    end
    Random.seed!(seed)
    p = randperm(len)
    ind = [p[r] for r in [b[i]:b[i+1] - 1 for i = 1:nfolds]]
    return [vector[i] for i in ind]
end

nfolds = 10
seed = 1
folds = [Vector{Int64}(undef,0) for _ in 1:nfolds]
## generate folds
for i in 1:length(fire_ids)
    println("starting fire: " * fire_ids[i])
    nf = gen_folds(training_boxes[training_boxes.fire_name .== fire_ids[i], :id], nfolds, seed)
    for j in 1:nfolds
        println("starting fold: " * string(j))
        f = []
        for k in 1:length(nf[j])
            xmn = training_boxes[training_boxes.id .== nf[j][k], :lower_left][1][1]
            ymn = training_boxes[training_boxes.id .== nf[j][k], :lower_left][1][2]
            xmx = xmn + bb
            ymx = ymn + bb
            f = vcat(f, findall(data.fire_name .== fire_ids[i] .&& data.x .>= xmn .&& data.x .< xmx .&&
                            data.y .>= ymn .&& data.y .< ymx))
        end
        folds[j] = vcat(folds[j], f)
    end
end

## reformat
folds_formatted = Tuple{Int, Int}[]
for i in 1:length(folds)
    folds_formatted = vcat(folds_formatted, (
        (folds[i,:][1]),
        (reduce(vcat, folds[1:end .!= i]))))
end

save_object("data/folds_formatted.jld2", folds_formatted)



##---------------------------------------------------------------
## For small scale metrics
##---------------------------------------------------------------

##---------------------------------------------------------------
## Fit and tune random forest model
##---------------------------------------------------------------

folds_formatted = load_object("data/folds_formatted.jld2");

data.holdout .= 1
data[vcat(folds_formatted[1][1], folds_formatted[1][2]), :holdout] .= 0

data.hs_class .= "not_high";
data[data.cbi .> 2.25, :hs_class] .= "high";

## target
target = Vector(data[:,:hs_class]);
target = coerce(target, OrderedFactor);
levels!(target, ["not_high", "high"]);

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


cols = [:clust_30_scaled, :mean_dens_30_scaled, :mean_area_30_scaled, :mean_ht_30_scaled, :em_scaled,
        :cwd_scaled, :slope_scaled, :tpi_scaled, :heat_load_scaled, :prev_sev_scaled, :avg_fuel_moisture_scaled,
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
## Tune model hyperparameters
##---------------------------------------------------------------

## kfold cross_validation
function kfold(model, fold_def, subsample_factor = 0.25, seed = 1)
    sum_ll = 0.0
    sum_auc = 0.0
    sum_mcc = 0.0
    sum_fscore = 0.0
    for i in 1:length(fold_def)
        println("Starting fold: " * string(i))
        Random.seed!(seed)
        sf1 = sample(folds_formatted[i][1], Int(round(length(folds_formatted[i][1]) * subsample_factor)))
        sf2 = sample(folds_formatted[i][2], Int(round(length(folds_formatted[i][2]) * subsample_factor)))
        MLJ.fit!(model, rows = sf2, force = true)
        yhat = MLJ.predict(model, rows = sf1)
        sum_ll += mean(log_loss(yhat, target[sf1]))
        sum_auc += auc(yhat, target[sf1])
        sum_mcc += mcc(mode.(yhat), target[sf1])
        sum_fscore += m(mode.(yhat), target[sf1])
    end
    return [sum_ll / length(fold_def), sum_auc / length(fold_def),
            sum_mcc / length(fold_def), sum_fscore / length(fold_def)]
end

## define f1score function
m = MulticlassFScore()

## define parameter range to test over
n_trees = [40]
max_depth = [20, 30, 40]
n_subfeatures = [1, 2, 3, 4]

p = collect(Base.product(n_trees, max_depth, n_subfeatures))
performance_data = DataFrame(log_loss = Vector{Float64}(undef, length(p)),
                             auc = Vector{Float64}(undef, length(p)),
                             mcc = Vector{Float64}(undef, length(p)),
                             fscore = Vector{Float64}(undef, length(p)));

## iterate through parameters and perform 10-fold CV
for i in 1:length(p)
    println("Starting iteration: " * string(i))
    forest_model = MLJDecisionTreeInterface.RandomForestClassifier(n_trees = p[i][1],
                                                                   max_depth = p[i][2],
                                                                   n_subfeatures = p[i][3])
    fm = machine(forest_model, features, target)
    performance_data[i,:] .= kfold(fm, folds_formatted, 0.25)
    CSV.write("data/rf_param_search.csv", performance_data)
end

## check results

CSV.write("data/model_tuning_30m.csv", performance_data)

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

save_object("data/forest_model.jld2", fm)

##---------------------------------------------------------------
## check out of sample prediction performance and compute feature importances
##---------------------------------------------------------------
fm = load_object("data/forest_model.jld2")

ho_list = findall(data.holdout .== 1);

yhat = MLJ.predict(fm, rows = ho_list);
m(mode.(yhat), data[ho_list, :hs_class])
mcc(mode.(yhat), data[ho_list, :hs_class])

auc(yhat, data[ho_list, :hs_class])

mean(log_loss(yhat, data[ho_list, :hs_class]))

confusion_matrix(mode.(yhat), data[ho_list, :hs_class])

##---------------------------------------------------------------
## Compute shapley metrics
##---------------------------------------------------------------
using Shapley: MonteCarlo

## create subset on which to compute shapley values
Random.seed!(1)
exp_list = sample(ho_list, 10000)
explain = copy(features[exp_list, :])

## compute shapley values
shap = shapley(η -> MLJ.predict(fm, η), MonteCarlo(CPUThreads(), 64), explain)
shap2 = shap
save_object("data/shapley.jld2", shap)

shap = load_object("data/shapley.jld2")

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
    xlims = [0, 0.05], frame = :box, legend = :none, color = :black, title = "Feature importance (high)",
    yflip = true)
savefig("plots/shapley_importance_high.pdf")

CSV.write("data/shapley_total_30m.csv", shapley_summary)

## summarize shapley values
shapley_summary = DataFrame(feature = names(features),
                            mean_not_high = Vector{Float64}(undef, ncol(features)),
                            sd_not_high = Vector{Float64}(undef, ncol(features)),
                            mean_high = Vector{Float64}(undef, ncol(features)),
                            sd_high = Vector{Float64}(undef, ncol(features)))

pi_sub = findall(data[exp_list, :own_type] .== "Private Industrial")
pf_sub = findall(data[exp_list, :own_type] .== "Federal")

length(shap)

for i in 1:length(shap)
    shp_mat = pdf(shap[i], ["not_high", "high"])
    shapley_summary[i, [2,4]] .= vec(mean(abs.(shp_mat[pi_sub]), dims = 1))
    shapley_summary[i, [3,5]] .= vec(std(abs.(shp_mat[pi_sub]), dims = 1))
end

sort(shapley_summary, order(:mean_high, rev = true))

data_plot = sort(shapley_summary, order(:mean_high, rev = true))
bar(data_plot.mean_high, orientation = :h, yticks = (1:16, data_plot.feature),
    xlims = [0, 0.07], frame = :box, legend = :none, color = :black, title = "Feature importance (high -- Private)",
    yflip = true)
savefig("plots/shapley_importance_high_private.pdf")

CSV.write("data/shapley_private_30m.csv", shapley_summary)

for i in 1:length(shap)
    shp_mat = pdf(shap[i], ["not_high","high"])
    shapley_summary[i, [2,4]] .= vec(mean(abs.(shp_mat[pf_sub]), dims = 1))
    shapley_summary[i, [3,5]] .= vec(std(abs.(shp_mat[pf_sub]), dims = 1))
end

sort(shapley_summary, order(:mean_high, rev = true))

data_plot = sort(shapley_summary, order(:mean_high, rev = true))
bar(data_plot.mean_high, orientation = :h, yticks = (1:16, data_plot.feature),
    xlims = [0, 0.05], frame = :box, legend = :none, color = :black, title = "Feature importance (high -- Public)",
    yflip = true)
savefig("plots/shapley_importance_high_public.pdf")

CSV.write("data/shapley_public_30m.csv", shapley_summary)

Random.seed!(2)
pf_ss = sample(pf_sub, length(pi_sub))

balanced_sample = vcat(pf_ss, pi_sub)

for i in 1:length(shap)
    shp_mat = pdf(shap[i], ["not_high",  "high"])
    shapley_summary[i, [2,4]] .= vec(mean(abs.(shp_mat[pf_ss]), dims = 1))
    shapley_summary[i, [3,5]] .= vec(std(abs.(shp_mat[pf_ss]), dims = 1))
end

sort(shapley_summary, order(:mean_high, rev = true))

data_plot = sort(shapley_summary, order(:mean_high, rev = true))
bar(data_plot.mean_high, orientation = :h, yticks = (1:16, data_plot.feature),
    xlims = [0, 0.05], frame = :box, legend = :none, color = :black, title = "Feature importance (high -- Balanced)",
    yflip = true)
savefig("plots/shapley_importance_high_balanced.pdf")

CSV.write("data/shapley_balanced_30m.csv", shapley_summary)




##---------------------------------------------------------------
## FOR LARGE SCALE METRICS
##
##
##
##---------------------------------------------------------------

##---------------------------------------------------------------
## Fit and tune random forest model
##---------------------------------------------------------------

folds_formatted = load_object("data/folds_formatted.jld2");

data.holdout .= 1
data[vcat(folds_formatted[1][1], folds_formatted[1][2]), :holdout] .= 0

data.hs_class .= "not_high";
data[data.cbi .> 2.25, :hs_class] .= "high";

## target
target = Vector(data[:,:hs_class]);
target = coerce(target, OrderedFactor);
levels!(target, ["not_high", "high"]);

## scale and center parameters

data.clust_180_scaled = (data.clust_180 .- mean(data.clust_180)) ./ std(data.clust_180)
data.mean_dens_180_scaled = (data.mean_dens_180 .- mean(data.mean_dens_180)) ./ std(data.mean_dens_180)
data.mean_area_180_scaled = (data.mean_area_180 .- mean(data.mean_area_180)) ./ std(data.mean_area_180)
data.mean_ht_180_scaled = (data.mean_ht_180 .- mean(data.mean_ht_180)) ./ std(data.mean_ht_180)

cols = [:clust_180_scaled, :mean_dens_180_scaled, :mean_area_180_scaled, :mean_ht_180_scaled, :em_scaled,
        :cwd_scaled, :slope_scaled, :tpi_scaled, :heat_load_scaled, :prev_sev_scaled, :avg_fuel_moisture_scaled,
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
## Tune model hyperparameters
##---------------------------------------------------------------

## define parameter range to test over
n_trees = [40]
max_depth = [20, 30, 40]
n_subfeatures = [1, 2, 3, 4]

p = collect(Base.product(n_trees, max_depth, n_subfeatures))
performance_data_180 = DataFrame(log_loss = Vector{Float64}(undef, length(p)),
                             auc = Vector{Float64}(undef, length(p)),
                             mcc = Vector{Float64}(undef, length(p)),
                             fscore = Vector{Float64}(undef, length(p)));

## iterate through parameters and perform 10-fold CV
for i in 1:length(p)
    println("Starting iteration: " * string(i))
    forest_model_180 = MLJDecisionTreeInterface.RandomForestClassifier(n_trees = p[i][1],
                                                                   max_depth = p[i][2],
                                                                   n_subfeatures = p[i][3])
    fm = machine(forest_model_180, features, target)
    performance_data_180[i,:] .= kfold(fm, folds_formatted, 0.25)
    CSV.write("data/rf_param_search_180.csv", performance_data_180)
end


## check results
performance_data_180

##---------------------------------------------------------------
## Fit model with tuned parameters
##---------------------------------------------------------------

forest_model_180 = MLJDecisionTreeInterface.RandomForestClassifier(n_trees = 40,
                                                               max_depth = 40,
                                                               n_subfeatures = 3)
fm_180 = machine(forest_model_180,
             features,
             target)

## fit model, takes awhile\
MLJ.fit!(fm_180, rows = vcat(folds_formatted[1][1], folds_formatted[1][2]))

save_object("data/forest_model_180.jld2", fm_180)

##---------------------------------------------------------------
## check out of sample prediction performance and compute feature importances
##---------------------------------------------------------------

fm_180 = load_object("data/forest_model_180.jld2")

ho_list = findall(data.holdout .== 1)

yhat_180 = MLJ.predict(fm_180, rows = ho_list)
m(mode.(yhat_180), data[ho_list, :hs_class])
mcc(mode.(yhat_180), data[ho_list, :hs_class])

auc(yhat_180, data[ho_list, :hs_class])

mean(log_loss(yhat_180, data[ho_list, :hs_class]))

confusion_matrix(mode.(yhat_180), data[ho_list, :hs_class])


##---------------------------------------------------------------
## Compute shapley metrics
##---------------------------------------------------------------

using Shapley: MonteCarlo

## create subset on which to compute shapley values
Random.seed!(2)
exp_list = sample(ho_list, 10000)
explain = copy(features[exp_list, :])

## compute shapley values
shap_180 = shapley(η -> MLJ.predict(fm_180, η), MonteCarlo(CPUThreads(), 64), explain)

save_object("data/shapley_180.jld2", shap_180)

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
    xlims = [0, 0.05], frame = :box, legend = :none, color = :black, title = "Feature importance (high)",
    yflip = true)
savefig("plots/shapley_importance_high_180.pdf")
CSV.write("data/shapley_total_180m.csv", shapley_summary_180)


## summarize shapley values
shapley_summary_180 = DataFrame(feature = names(features),
                                mean_not_high = Vector{Float64}(undef, ncol(features)),
                                sd_not_high = Vector{Float64}(undef, ncol(features)),
                                mean_high = Vector{Float64}(undef, ncol(features)),
                                sd_high = Vector{Float64}(undef, ncol(features)))

pi_sub = findall(data[exp_list, :own_type] .== "Private Industrial")
pf_sub = findall(data[exp_list, :own_type] .== "Federal")

for i in 1:length(shap_180)
    shp_mat = pdf(shap_180[i], ["not_high", "high"])
    shapley_summary_180[i, [2,4]] .= vec(mean(abs.(shp_mat[pi_sub]), dims = 1))
    shapley_summary_180[i, [3,5]] .= vec(std(abs.(shp_mat[pi_sub]), dims = 1))
end

sort(shapley_summary_180, order(:mean_high, rev = true))

data_plot = sort(shapley_summary_180, order(:mean_high, rev = true))
bar(data_plot.mean_high, orientation = :h, yticks = (1:16, data_plot.feature),
    xlims = [0, 0.07], frame = :box, legend = :none, color = :black, title = "Feature importance (high -- Private)",
    yflip = true)
savefig("plots/shapley_importance_high_private_180.pdf")

CSV.write("data/shapley_private_180m.csv", shapley_summary_180)


for i in 1:length(shap_180)
    shp_mat = pdf(shap_180[i], ["not_high","high"])
    shapley_summary_180[i, [2,4]] .= vec(mean(abs.(shp_mat[pf_sub]), dims = 1))
    shapley_summary_180[i, [3,5]] .= vec(std(abs.(shp_mat[pf_sub]), dims = 1))
end

sort(shapley_summary_180, order(:mean_high, rev = true))

data_plot = sort(shapley_summary_180, order(:mean_high, rev = true))
bar(data_plot.mean_high, orientation = :h, yticks = (1:16, data_plot.feature),
    xlims = [0, 0.05], frame = :box, legend = :none, color = :black, title = "Feature importance (high -- Public)",
    yflip = true)
savefig("plots/shapley_importance_high_public_180.pdf")
CSV.write("data/shapley_public_180m.csv", shapley_summary_180)

Random.seed!(2)
pf_ss = sample(pf_sub, length(pi_sub))

balanced_sample = vcat(pf_ss, pi_sub)

for i in 1:length(shap_180)
    shp_mat = pdf(shap_180[i], ["not_high",  "high"])
    shapley_summary_180[i, [2,4]] .= vec(mean(abs.(shp_mat[pf_ss]), dims = 1))
    shapley_summary_180[i, [3,5]] .= vec(std(abs.(shp_mat[pf_ss]), dims = 1))
end

sort(shapley_summary_180, order(:mean_high, rev = true))

data_plot = sort(shapley_summary_180, order(:mean_high, rev = true))
bar(data_plot.mean_high, orientation = :h, yticks = (1:16, data_plot.feature),
    xlims = [0, 0.05], frame = :box, legend = :none, color = :black, title = "Feature importance (high -- Balanced)",
    yflip = true)
savefig("plots/shapley_importance_high_balanced.pdf")
CSV.write("data/shapley_balanced_180m.csv", shapley_summary_180)


##---------------------------------------------------------------
## Partial dependency plots
##---------------------------------------------------------------

function pdp(;var = [:mean_dens_30], npts = 100, sp = 0.5, st = -1, var2_levels = [0], pts = true)

    bdf = DataFrame(fire_name = ["sugar", "walker", "north_complex", "sheep"],
                    clust_30 = zeros(4),
                    mean_dens_30 = zeros(4),
                    mean_area_30 = zeros(4),
                    sd_area_30 = zeros(4),
                    mean_ht_30 = zeros(4),
                    mean_frac_30 = zeros(4),
                    em = zeros(4),
                    cwd = zeros(4),
                    slope = zeros(4),
                    tpi = zeros(4),
                    heat_load = zeros(4),
                    prev_sev = zeros(4),
                    avg_fuel_moisture = zeros(4),
                    hdw = zeros(4),
                    avg_fuel_temp = zeros(4));
    varlist = [:fire_name, :clust_30_scaled, :mean_dens_30_scaled, :mean_area_30_scaled,
               :mean_ht_30_scaled, :em_scaled, :cwd_scaled, :slope_scaled,
               :tpi_scaled, :heat_load_scaled, :prev_sev_scaled, :avg_fuel_moisture_scaled, :hdw_scaled]

    if st == -1
        st = round((maximum(data[:,var[1]]) - minimum(data[:,var[1]])) / npts, digits = 2)
    end

    if length(var) > 1
        pd = DataFrame(Base.product(collect(range(0, round(maximum(data.mean_dens_30)), length = 100)),
                            var2_levels));
        rename!(pd, var)
        for v in varlist[2:length(varlist)][varlist[2:length(varlist)] .∉ [var]]
            pd[:,v] .= median(data[:,v])
        end
    else
        pd = DataFrame()
        pd[:,var[1]] = collect(range(round(minimum(data[:,var[1]])), round(maximum(data[:,var[1]])), length = npts))
        for v in varlist[2:length(varlist)][varlist[2:length(varlist)] .!= var[1]]
            pd[:,v] .= median(data[:,v])
        end
    end

    ## coerce to correct scitype
    for i in 1:ncol(pd)
        pd[:,i]  = coerce(float.(identity.(pd[:,i])), Continuous)
    end

    pd.fire_name .= "dixie";
    pd = vcat(pd, bdf);

    pd.fire_name = coerce(pd.fire_name, OrderedFactor);
    pd = pd[:, [:fire_name, :clust_30_scaled, :mean_dens_30_scaled, :mean_area_30_scaled,
               :mean_ht_30_scaled, :em_scaled, :cwd_scaled, :slope_scaled,
               :tpi_scaled, :heat_load_scaled, :prev_sev_scaled, :avg_fuel_moisture_scaled, :hdw_scaled]];

    yt = MLJ.predict(fm, pd);

    prs = pdf(yt, ["not_high", "high"]);

    # layout = @layout [a
    #                   b{1.0w,0.8h}]

    if length(var) > 1

        colors = ["#ffeda0", "#fec44f", "#f03b20"]

        np = npts * length(var2_levels)

        phigh = Plots.plot(layout = layout, link = :both, size = (500, 500), margin = -10Plots.px)
        if pts
            plot!(pd[1:np, var[1]], prs[1:np,2], subplot = 2, group = pd[1:np, var[2]],
                  seriestype = :scatter, markercolor = :match, color_palette = colors)
        end
        dx = maximum(data[:, var[1]]) / 150
        binedges = 0:dx:maximum(data[:, var[1]])
        spr = copy(data[data.own_type .== "Private Industrial", var[1]])
        sf = copy(data[data.own_type .== "Federal", var[1]])
        aw = histcounts(spr[sample(1:length(spr), 500000)], binedges)
        bw = histcounts(sf[sample(1:length(sf), 500000)], binedges)
        for i in 1:length(unique(pd[:,var[2]]))-1
            lm = loess(pd[pd[!,var[2]] .== unique(pd[!,var[2]])[i], var[1]], prs[pd[!,var[2]] .== unique(pd[!,var[2]])[i],2], span = 0.5)
            nx = range(extrema(pd[pd[!,var[2]] .== unique(pd[!,var[2]])[i], var[1]])...; step = 5)
            ny = Loess.predict(lm, nx)
            plot!(nx, ny, subplot = 2, color = colors[i], linewidth = 3, legend = :none, frame = :box)
        end

        xlabel!(string(var[1]))
        bar!(binedges, aw, subplot = 1, bar_width = dx, alpha = 0.3,
             frame = :none, color = :blue, legend = :none, linewidth = 0)
        bar!(binedges, bw, subplot = 1, bar_width = dx, alpha = 0.3,
             frame = :none, color = :green, legend = :none, linewidth = 0, title = "High")

        # plow = Plots.plot(layout = layout, link = :both, size = (500, 500), margin = -10Plots.px)
        # if pts
        #     plot!(pd[1:np, var[1]], prs[1:np,1], group = pd[1:np, var[2]], subplot = 2,
        #           seriestype = :scatter, markercolor = :match, color_palette = colors)
        # end
        # for i in 1:length(unique(pd[:,var[2]]))-1
        #     lm = loess(pd[pd[!,var[2]] .== unique(pd[!,var[2]])[i], var[1]], prs[pd[!,var[2]] .== unique(pd[!,var[2]])[i],1], span = 0.5)
        #     nx = range(extrema(pd[pd[!,var[2]] .== unique(pd[!,var[2]])[i], var[1]])...; step = 5)
        #     ny = Loess.predict(lm, nx)
        #     plot!(nx, ny, subplot = 2, color = colors[i], linewidth = 3, legend = :none, frame = :box)
        # end
        # xlabel!(string(var[1]))
        # bar!(binedges, aw, subplot = 1, bar_width = dx, alpha = 0.3,
        #      frame = :none, color = :blue, legend = :none, linewidth = 0)
        # bar!(binedges, bw, subplot = 1, bar_width = dx, alpha = 0.3,
        #      frame = :none, color = :green, legend = :none, linewidth = 0, title = "Low")

    else

        dx = maximum(data[:, var[1]]) / 150
        binedges = 0:dx:maximum(data[:, var[1]])
        spr = copy(data[data.own_type .== "Private Industrial", var[1]])
        sf = copy(data[data.own_type .== "Federal", var[1]])
        aw = histcounts(spr[sample(1:length(spr), 500000)], binedges)
        bw = histcounts(sf[sample(1:length(sf), 500000)], binedges)

        phigh = Plots.plot(layout = layout, link = :both, size = (400, 400), margin = -10Plots.px)
        if pts
            plot!(pd[1:npts, var[1]], prs[1:npts,2], seriestype = :scatter, color = :black, subplot = 2)
        end
        lm = loess(pd[1:npts, var[1]], prs[1:npts,2], span = sp)
        nx = range(extrema(pd[1:npts, var[1]])...; step = st)
        ny = Loess.predict(lm, nx)
        plot!(nx, ny, color = :red, linewidth = 3, legend = :none, frame = :box,
              ylims = [round(minimum(prs), digits = 2), round(maximum(prs), digits = 2)], subplot = 2)
        xlabel!(string(var[1]))
        bar!(binedges, aw, subplot = 1, bar_width = dx, alpha = 0.3,
             frame = :none, color = :blue, legend = :none, linewidth = 0)
        bar!(binedges, bw, subplot = 1, bar_width = dx, alpha = 0.3,
             frame = :none, color = :green, legend = :none, linewidth = 0, title = "High")

        # plow = Plots.plot(layout = layout, link = :both, size = (400, 400), margin = -10Plots.px)
        # if pts
        #     plot!(pd[1:npts, var[1]], prs[1:npts,1], seriestype = :scatter, color = :black, subplot = 2)
        # end
        # lm = loess(pd[1:npts, var[1]], prs[1:npts,1], span = sp)
        # nx = range(extrema(pd[1:npts, var[1]])...; step = st)
        # ny = Loess.predict(lm, nx)
        # plot!(nx, ny, color = :green, linewidth = 3, legend = :none, frame = :box,
        #       ylims = [round(minimum(prs), digits = 2), round(maximum(prs), digits = 2)], subplot = 2)
        # xlabel!(string(var[1]))
        # bar!(binedges, aw, subplot = 1, bar_width = dx, alpha = 0.3,
        #      frame = :none, color = :blue, legend = :none, linewidth = 0)
        # bar!(binedges, bw, subplot = 1, bar_width = dx, alpha = 0.3,
        #      frame = :none, color = :green, legend = :none, linewidth = 0, title = "Low")

    end

    return Plots.plot(phigh)

end

## forest structure
pdp(var = [:mean_dens_30], pts = false, npts = 500)
savefig("../plots/pdp_mean_dens.pdf")
pdp(var = [:clust_30], pts = false, npts = 500)
savefig("../plots/pdp_clust.pdf")
pdp(var = [:mean_area_30], pts = false)
savefig("../plots/pdp_mean_area.pdf")
pdp(var = [:sd_area_30], pts = false)
savefig("../plots/pdp_sd_area.pdf")
pdp(var = [:em], sp = 0.3, pts = false)
savefig("../plots/pdp_em.pdf")
pdp(var = [:mean_frac_30], sp = 0.3, pts = false)

## weather
pdp(var = [:hdw], pts = false, npts = 500)
savefig("../plots/pdp_hdw.pdf")
pdp(var = [:avg_fuel_moisture], sp = 0.3, pts = false)
savefig("../plots/pdp_fuel_moisture.pdf")
pdp(var = [:avg_fuel_temp], sp = 0.3, pts = false)
savefig("../plots/pdp_avg_fuel_temp.pdf")
pdp(var = [:prev_sev], pts = false)
savefig("../plots/pdp_prev_sev.pdf")

## topography
pdp(var = [:cwd], pts = false)
savefig("../plots/pdp_cwd.pdf")
pdp(var = [:slope], pts = false)
savefig("../plots/pdp_slope.pdf")
pdp(var = [:tpi], pts = false)
savefig("../plots/pdp_tpi.pdf")
pdp(var = [:heat_load], pts = false)
savefig("../plots/pdp_heat_load.pdf")

pdp(var = [:mean_dens_30, :hdw], var2_levels = [5, 75, 150], pts = false)
savefig("../plots/pdp_dens_hdw.pdf")

pdp(var = [:mean_dens_30, :prev_sev], var2_levels = [0, 1, 2], pts = false)


pdp(var = [:mean_dens_30, :em], var2_levels = [10, 15, 20], pts = false)

pdp(var = [:mean_dens, :prev_sev], var2_levels = [0.5, 1.5, 2.5], pts = false)
savefig("../plots/pdp_dens_prevsev.pdf")

pdp(var = [:hdw, :mean_dens_30], var2_levels = [50, 150, 250], pts = false)
savefig("../plots/pdp_hdw_dens.pdf")

##---------------------------------------------------------------
## predictions when only varying forest structure
##---------------------------------------------------------------

fake_features = copy(features);
fake_features.cwd .= mean(features.cwd);
fake_features.slope .= mean(features.slope);
fake_features.tpi .= mean(features.tpi);
fake_features.heat_load .= mean(features.heat_load);
fake_features.prev_sev .= mean(features.prev_sev);
fake_features.avg_fuel_moisture .= mean(features.avg_fuel_moisture);
fake_features.hdw .= mean(features.hdw);
fake_features.avg_fuel_temp .= mean(features.avg_fuel_temp);

fake_predict = MLJ.predict(fm, fake_features[ho_list, :]);

probs = pdf(yhat, ["low-none", "med", "high"]);

pi_list = findall(data[ho_list, :own_type] .== "Private Industrial" .&& data[ho_list, :cwd] .< 1500)
pf_list = findall(data[ho_list, :own_type] .== "Federal" .&& data[ho_list, :cwd] .< 1500)

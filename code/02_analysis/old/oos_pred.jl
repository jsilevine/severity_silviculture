
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

testing_boxes = boxes_full[boxes_full.fire_name .== "sheep" .|| boxes_full.fire_name .== "walker",:]
training_boxes = boxes_full[boxes_full.fire_name .!= "sheep" .&& boxes_full.fire_name .!= "walker",:]

##---------------------------------------------------------------
## Generate folds for CV
##---------------------------------------------------------------

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

fids = unique(training_boxes.fire_name)
## generate folds
for i in 1:length(fids)
    println("starting fire: " * fids[i])
    nf = gen_folds(training_boxes[training_boxes.fire_name .== fids[i], :id], nfolds, seed)
    for j in 1:nfolds
        println("starting fold: " * string(j))
        f = []
        for k in 1:length(nf[j])
            xmn = training_boxes[training_boxes.id .== nf[j][k], :lower_left][1][1]
            ymn = training_boxes[training_boxes.id .== nf[j][k], :lower_left][1][2]
            xmx = xmn + bb
            ymx = ymn + bb
            f = vcat(f, findall(data.fire_name .== fids[i] .&& data.x .>= xmn .&& data.x .< xmx .&&
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

save_object("data/folds_formatted_oos.jld2", folds_formatted)

minimum(data.cwd)
histogram(data.cwd)
quantile(data.cwd, 0.05)
quantile(data.cwd, 0.95)
##---------------------------------------------------------------
## FOR LARGE SCALE METRICS
##
##
##---------------------------------------------------------------

##---------------------------------------------------------------
## Fit and tune random forest model
##---------------------------------------------------------------

folds_formatted = load_object("data/folds_formatted_oos.jld2");

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
data.em_scaled = (data.em .- mean(data.em)) ./ std(data.em);
data.cwd_scaled = (data.cwd .- mean(data.cwd)) ./ std(data.cwd);
data.slope_scaled = (data.slope .- mean(data.slope)) ./ std(data.slope);
data.tpi_scaled = (data.tpi .- mean(data.tpi)) ./ std(data.tpi);
data.heat_load_scaled = (data.heat_load .- mean(data.heat_load)) ./ std(data.heat_load);
data.prev_sev_scaled = (data.prev_sev .- mean(filter(!isnan, data.prev_sev))) ./ std(filter(!isnan, data.prev_sev));
data.avg_fuel_moisture_scaled = (data.avg_fuel_moisture .- mean(data.avg_fuel_moisture)) ./ std(data.avg_fuel_moisture);
data.hdw_scaled = (data.hdw .- mean(data.hdw)) ./ std(data.hdw);

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

save_object("data/forest_model_180_oos.jld2", fm_180)

##---------------------------------------------------------------
## check out of sample prediction performance and compute feature importances
##---------------------------------------------------------------

fm_180 = load_object("data/forest_model_180_oos.jld2")

ho_list = findall(data.holdout .== 1)

yhat_180 = MLJ.predict(fm_180, rows = ho_list)
m(mode.(yhat_180), data[ho_list, :hs_class])
mcc(mode.(yhat_180), data[ho_list, :hs_class])

auc(yhat_180, data[ho_list, :hs_class])

mean(log_loss(yhat_180, data[ho_list, :hs_class]))

confusion_matrix(mode.(yhat_180), data[ho_list, :hs_class])

##---------------------------------------------------------------
## Removing CWD
##---------------------------------------------------------------

select!(features, Not(:cwd_scaled))

fm_180_nocwd = machine(forest_model_180,
                       features,
                       target)

## fit model, takes awhile\
MLJ.fit!(fm_180_nocwd, rows = vcat(folds_formatted[1][1], folds_formatted[1][2]))

ho_list = findall(data.holdout .== 1)

yhat_180_nocwd = MLJ.predict(fm_180_nocwd, rows = ho_list)

auc(yhat_180_nocwd, data[ho_list, :hs_class])

confusion_matrix(mode.(yhat_180_nocwd), data[ho_list, :hs_class])

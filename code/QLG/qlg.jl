using Pkg




using Base: index_shape
using CSV, DataFrames, MLJ, MLJDecisionTreeInterface,
    DecisionTree, Plots, StatsBase, Random, MLJBase,
    JLD2, Shapley, Shapley.ComputationalResources,
    NaNStatistics, Loess

## load in data
data = CSV.read("../../data/complete_data.csv", DataFrame);
data = data[completecases(data),:];

##---------------------------------------------------------------
## create grid for kfold cross validation
##---------------------------------------------------------------

bb = CSV.read("../data/autocor_scale.csv", DataFrame)[1,2];

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
CSV.write("../data/boxes_for_kfold.csv", boxes_full)


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

CSV.write("../data/testing_boxes.csv", testing_boxes)
CSV.write("../data/training_boxes.csv", training_boxes)


##---------------------------------------------------------------
## Generate folds for CV
##---------------------------------------------------------------
training_boxes = CSV.read("../data/training_boxes.csv", DataFrame)
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

save_object("../data/folds_formatted.jld2", folds_formatted)

##---------------------------------------------------------------
## Fit and tune random forest model
##---------------------------------------------------------------

folds_formatted = load_object("../../data/folds_formatted.jld2");

data.holdout .= 1
data[vcat(folds_formatted[1][1], folds_formatted[1][2]), :holdout] .= 0

data.hs_class .= "not-high";
data[data.CBI_bc .> 2.25, :hs_class] .= "high";

## target
target = Vector(data[:,:hs_class]);
target = coerce(target, OrderedFactor);
levels!(target, ["not-high", "high"]);

cols = [:clust, :mean_dens, :mean_area, :sd_area, :em,
        :elevation, :slope, :tpi, :heat_load, :prev_sev, :avg_fuel_moisture,
        :hdw, :avg_fuel_temp, :percent_open];
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
    CSV.write("../../data/rf_param_search.csv", performance_data)
end

## check results
performance_data

##---------------------------------------------------------------
## Fit model with tuned parameters
##---------------------------------------------------------------

forest_model = MLJDecisionTreeInterface.RandomForestClassifier(n_trees = 40,
                                                               max_depth = 40,
                                                               n_subfeatures = 4)
fm = machine(forest_model,
             features,
             target)

## fit model, takes awhile
MLJ.fit!(fm, rows = vcat(folds_formatted[1][1], folds_formatted[1][2]))

save_object("../../data/forest_model_qlg.jld2", fm)

##---------------------------------------------------------------
## check out of sample prediction performance and compute feature importances
##---------------------------------------------------------------
fm = load_object("../../data/forest_model_qlg.jld2")

ho_list = findall(data.holdout .== 1)

yhat = MLJ.predict(fm, rows = ho_list)

yhat_mode = MLJ.predict_mode(fm, rows = ho_list)
yhat_mode = categorical(yhat_mode)
m(mode.(yhat), data[ho_list, :hs_class])
mcc(mode.(yhat), data[ho_list, :hs_class])

auc(yhat, data[ho_list, :hs_class])
mean(log_loss(yhat, data[ho_list, :hs_class]))

confusion_matrix(mode.(yhat), data[ho_list, :hs_class])


##---------------------------------------------------------------
## Predictions from new data
##---------------------------------------------------------------

## get proposed data
proposed_dfpz_west = CSV.read("../../data/qlg/proposed_txs_forest_structure_unmodified/proposedDfpz_forestStruc_west.csv", DataFrame)
proposed_dfpz_transition = CSV.read("../../data/qlg/proposed_txs_forest_structure_unmodified/proposedDfpz_forestStruc_transition.csv", DataFrame)
proposed_dfpz_east = CSV.read("../../data/qlg/proposed_txs_forest_structure_unmodified/proposedDfpz_forestStruc_east.csv", DataFrame)

proposed_dfpz_west.zone .= "west"
proposed_dfpz_transition.zone .= "transition"
proposed_dfpz_east.zone .= "east"

proposed_dfpz = vcat(proposed_dfpz_west[:,[:x_orig, :y_orig, :zone]],
                     proposed_dfpz_transition[:,[:x_orig, :y_orig, :zone]],
                     proposed_dfpz_east[:,[:x_orig, :y_orig, :zone]])

proposed_matrix_west = CSV.read("../../data/qlg/proposed_txs_forest_structure_unmodified/proposedMatrix_forestStruc_west.csv", DataFrame)
proposed_matrix_transition = CSV.read("../../data/qlg/proposed_txs_forest_structure_unmodified/proposedMatrix_forestStruc_transition.csv", DataFrame)
proposed_matrix_east = CSV.read("../../data/qlg/proposed_txs_forest_structure_unmodified/proposedMatrix_forestStruc_east.csv", DataFrame)

proposed_matrix_west.zone .= "west"
proposed_matrix_transition.zone .= "transition"
proposed_matrix_east.zone .= "east"

proposed_matrix = vcat(proposed_matrix_west[:,[:x_orig, :y_orig, :zone]],
                     proposed_matrix_transition[:,[:x_orig, :y_orig, :zone]],
                     proposed_matrix_east[:,[:x_orig, :y_orig, :zone]])

proposed_dfpz.qlg_treatment_type .= "dfpz"
proposed_matrix.qlg_treatment_type .= "matrix"

proposed_qlg_treatments = vcat(proposed_dfpz, proposed_matrix)

rename!(proposed_qlg_treatments, :x_orig => :x, :y_orig => :y)

full_data_qlg = leftjoin(data, proposed_qlg_treatments, on = [:x, :y])

full_data_qlg[ismissing.(full_data_qlg.qlg_treatment_type), :qlg_treatment_type] .= "none"
full_data_qlg[ismissing.(full_data_qlg.zone), :zone] .= "none"

forest_structure_west = CSV.read("../../data/qlg/existing_txs_forest_structure/qlg_forestStructureTreated_west_erasePre18Fires.csv", DataFrame)
forest_structure_transition = CSV.read("../../data/qlg/existing_txs_forest_structure/qlg_forestStructureTreated_transition_erasePre18Fires.csv", DataFrame)
forest_structure_east = CSV.read("../../data/qlg/existing_txs_forest_structure/qlg_forestStructureTreated_east_erasePre18Fires.csv", DataFrame)

rename!(forest_structure_west, :clustering_index => :clust, :mean_stem_density_per_ha => :mean_dens,
        :mean_gap_area_m => :mean_area, :sd_gap_area_m => :sd_area, :ladder_fuels_metric => :em)
rename!(forest_structure_transition, :clustering_index => :clust, :mean_stem_density_per_ha => :mean_dens,
        :mean_gap_area_m => :mean_area, :sd_gap_area_m => :sd_area, :ladder_fuels_metric => :em)
rename!(forest_structure_east, :clustering_index => :clust, :mean_stem_density_per_ha => :mean_dens,
        :mean_gap_area_m => :mean_area, :sd_gap_area_m => :sd_area, :ladder_fuels_metric => :em)

cnames = [:clust, :mean_dens, :mean_area, :sd_area, :em, :percent_open]
forest_structure_west = forest_structure_west[completecases(forest_structure_west[:, cnames]), cnames]
forest_structure_transition = forest_structure_transition[completecases(forest_structure_transition[:, cnames]), cnames]
forest_structure_east = forest_structure_east[completecases(forest_structure_east[:, cnames]), cnames]


cols = [:clust, :mean_dens, :mean_area, :sd_area, :em,
        :elevation, :slope, :tpi, :heat_load, :prev_sev, :avg_fuel_moisture,
        :hdw, :avg_fuel_temp, :percent_open];
features = DataFrame(fire_name = coerce(full_data_qlg[:,:fire_name], OrderedFactor));

## coerce to correct scitype
for i in 1:length(cols)
    features[:,cols[i]]  = coerce(float.(identity.(full_data_qlg[:,cols[i]])), Continuous)
end

ynew_baseline = MLJ.predict(fm, features)
ynew_baseline_mode = mode.(ynew_baseline)
prob_high_baseline = pdf(ynew_baseline, ["low-none", "med", "high"])[:,3]


niter = 20
results = Array{Any}(missing, nrow(full_data_qlg), niter)


##---------------------------------------------------------------
## Get baseline for untreated areas
##---------------------------------------------------------------

temp_data = copy(full_data_qlg)
complete_untrt_indx = completecases(temp_data) .&& temp_data.qlg_treatment_type .== "none"

temp_data = temp_data[complete_untrt_indx,:]
cols = [:clust, :mean_dens, :mean_area, :sd_area, :em,
        :elevation, :slope, :tpi, :heat_load, :prev_sev, :avg_fuel_moisture,
        :hdw, :avg_fuel_temp, :percent_open];
features = DataFrame(fire_name = coerce(temp_data[:,:fire_name], OrderedFactor));

## coerce to correct scitype
for i in 1:length(cols)
    println(string(cols[i]))
    features[:,cols[i]]  = coerce(float.(identity.(temp_data[:,cols[i]])), Continuous)
end

p = MLJ.predict(fm, features)

p_pdf = pdf(p, ["not-high", "high"])[:,2]
p_mode = mode.(p)
for i in 1:size(results)[2]
    println(i)
    results[complete_untrt_indx,i] .= p_pdf
    results_mode[complete_untrt_indx,i] .= p_mode
end


baseline = DataFrame(x = full_data_qlg.x,
                     y = full_data_qlg.y)
baseline.pr .= 0.0
baseline[complete_untrt_indx, :pr] .= p

p = Nothing

##---------------------------------------------------------------
## get baseline for treated areas
##---------------------------------------------------------------
temp_data = copy(full_data_qlg)
complete_trt_indx = completecases(temp_data) .&& temp_data.qlg_treatment_type .!= "none"
res = full_data_qlg[full_data_qlg.fire_name .== "sheep", :]
res = res[1:2,:] ## add one datapoint from sheep
temp_data = temp_data[complete_trt_indx,:]
temp_data = vcat(temp_data, res)

cols = [:clust, :mean_dens, :mean_area, :sd_area, :em,
        :elevation, :slope, :tpi, :heat_load, :prev_sev, :avg_fuel_moisture,
        :hdw, :avg_fuel_temp, :percent_open];
features = DataFrame(fire_name = coerce(temp_data[:,:fire_name], OrderedFactor));

## coerce to correct scitype
for i in 1:length(cols)
    println(string(cols[i]))
    features[:,cols[i]]  = coerce(float.(identity.(temp_data[:,cols[i]])), Continuous)
end

baseline_trt = MLJ.predict(fm, features)
baseline_trt = baseline_trt[1:length(baseline_trt)-2]
baseline[]

baseline[complete_trt_indx, :pr] .= baseline_treated[:,:pr]
baseline.qlg_treatment_type .= "none"
baseline[complete_trt_indx, :qlg_treatment_type] .= "treated"

mean(baseline[baseline.qlg_treatment_type .== "none", :pr])
mean(baseline[baseline.qlg_treatment_type .== "treated", :pr])
CSV.write("../../data/qlg/baseline_treated.csv", baseline)

CSV.write("../../data/qlg/complete_data_qlg.csv", full_data_qlg)

##---------------------------------------------------------------
## do predictions on fake data
##---------------------------------------------------------------


full_data_qlg = CSV.read("../../data/qlg/complete_data_qlg.csv", DataFrame)


for i in 2:niter

    println("starting iteration: " * string(i))
    temp_data = copy(full_data_qlg)

    cnames = [:clust, :mean_dens, :mean_area, :sd_area, :em, :percent_open]

    ## for west
    indx = [1:1:nrow(forest_structure_west);][forest_structure_west.percent_open .<=
        quantile(forest_structure_west.percent_open, 0.75)]
    temp_data[temp_data.zone .== "west" .&& temp_data.qlg_treatment_type .== "dfpz", cnames] .=
        forest_structure_west[sample(1:nrow(forest_structure_west),
                                     nrow(temp_data[temp_data.zone .== "west" .&&
                                         temp_data.qlg_treatment_type .== "dfpz", :])), cnames]

    for k in 1:3
        list = []
        for j in findall(temp_data.zone .== "west" .&& temp_data.qlg_treatment_type .== "dfpz")
            if temp_data[j, :mean_dens] > full_data_qlg[j, :mean_dens]
                push!(list, j)
            end
        end
        temp_data[list, cnames] .= forest_structure_west[sample(1:nrow(forest_structure_west),
                                                                length(list)), cnames]
    end

    temp_data[temp_data.zone .== "west" .&& temp_data.qlg_treatment_type .== "matrix", cnames] .=
        forest_structure_west[sample(indx,
                                     nrow(temp_data[temp_data.zone .== "west" .&&
                                         temp_data.qlg_treatment_type .== "matrix", :])), cnames]
    for k in 1:3
        list = []
        for j in findall(temp_data.zone .== "west" .&& temp_data.qlg_treatment_type .== "matrix")
            if temp_data[j, :mean_dens] > full_data_qlg[j, :mean_dens]
                push!(list, j)
            end
        end
        temp_data[list, cnames] .= forest_structure_west[sample(indx, length(list)), cnames]
    end


    ## for transition
    indx = [1:1:nrow(forest_structure_transition);][forest_structure_transition.percent_open .<=
        quantile(forest_structure_transition.percent_open, 0.75)]
    temp_data[temp_data.zone .== "transition" .&& temp_data.qlg_treatment_type .== "dfpz", cnames] .=
        forest_structure_transition[sample(1:nrow(forest_structure_transition),
                                     nrow(temp_data[temp_data.zone .== "transition" .&&
                                         temp_data.qlg_treatment_type .== "dfpz", :])), cnames]

    for k in 1:3
        list = []
        for j in findall(temp_data.zone .== "transition" .&& temp_data.qlg_treatment_type .== "dfpz")
            if temp_data[j, :mean_dens] > full_data_qlg[j, :mean_dens]
                push!(list, j)
            end
        end
        temp_data[list, cnames] .= forest_structure_transition[sample(1:nrow(forest_structure_transition),
                                                                      length(list)), cnames]
    end

    temp_data[temp_data.zone .== "transition" .&& temp_data.qlg_treatment_type .== "matrix", cnames] .=
        forest_structure_transition[sample(indx,
                                     nrow(temp_data[temp_data.zone .== "transition" .&&
                                         temp_data.qlg_treatment_type .== "matrix", :])), cnames]

    for k in 1:3
        list = []
        for j in findall(temp_data.zone .== "transition" .&& temp_data.qlg_treatment_type .== "matrix")
            if temp_data[j, :mean_dens] > full_data_qlg[j, :mean_dens]
                push!(list, j)
            end
        end
        temp_data[list, cnames] .= forest_structure_transition[sample(indx,
                                                                      length(list)), cnames]
    end

    ## for east
    indx = [1:1:nrow(forest_structure_east);][forest_structure_east.percent_open .<=
        quantile(forest_structure_east.percent_open, 0.75)]
    temp_data[temp_data.zone .== "east" .&& temp_data.qlg_treatment_type .== "dfpz", cnames] .=
        forest_structure_east[sample(1:nrow(forest_structure_east),
                                     nrow(temp_data[temp_data.zone .== "east" .&&
                                         temp_data.qlg_treatment_type .== "dfpz", :])), cnames]

    for k in 1:3
        list = []
        for j in findall(temp_data.zone .== "east" .&& temp_data.qlg_treatment_type .== "dfpz")
            if temp_data[j, :mean_dens] > full_data_qlg[j, :mean_dens]
                push!(list, j)
            end
        end
        temp_data[list, cnames] .= forest_structure_east[sample(1:nrow(forest_structure_east),
                                                                length(list)), cnames]
    end

    temp_data[temp_data.zone .== "east" .&& temp_data.qlg_treatment_type .== "matrix", cnames] .=
        forest_structure_east[sample(indx,
                                     nrow(temp_data[temp_data.zone .== "east" .&&
                                         temp_data.qlg_treatment_type .== "matrix", :])), cnames]
    for k in 1:3
        list = []
        for j in findall(temp_data.zone .== "east" .&& temp_data.qlg_treatment_type .== "matrix")
            if temp_data[j, :mean_dens] > full_data_qlg[j, :mean_dens]
                push!(list, j)
            end
        end
        temp_data[list, cnames] .= forest_structure_east[sample(indx, length(list)), cnames]
    end

    ## if density increased due to treatment, revert to original value
    for j in findall(temp_data.qlg_treatment_type .!= "none")
        if temp_data[j, :mean_dens] > full_data_qlg[j, :mean_dens]
            temp_data[j, :] = full_data_qlg[j, :]
        end
    end

    complete_trt_indx = completecases(temp_data) .&& temp_data.qlg_treatment_type .!= "none"
    res = full_data_qlg[full_data_qlg.fire_name .== "sheep", :]
    res = res[1:2,:] ## add one datapoint from sheep
    temp_data = temp_data[complete_trt_indx,:]
    temp_data = vcat(temp_data, res)

    cols = [:clust, :mean_dens, :mean_area, :sd_area, :em,
            :elevation, :slope, :tpi, :heat_load, :prev_sev, :avg_fuel_moisture,
            :hdw, :avg_fuel_temp, :percent_open];
    features = DataFrame(fire_name = coerce(temp_data[:,:fire_name], OrderedFactor));

    ## coerce to correct scitype
    for i in 1:length(cols)
        println(string(cols[i]))
        features[:,cols[i]]  = coerce(float.(identity.(temp_data[:,cols[i]])), Continuous)
    end

    println("starting predictions")
    p = MLJ.predict(fm, features)
    p = p[1:length(p)-2]

    results[complete_trt_indx,i] .= pdf(p, ["not-high", "high"])[:,2]

    CSV.write("../../data/qlg/results_prob.csv", DataFrame(results, :auto))

end


ynew_mode = mode.(ynew)
prob_high = pdf(ynew, ["low-none", "med", "high"])[:,3]

out = DataFrame(x = full_data_qlg.x, y = full_data_qlg.y,
                sev_pred = ynew_mode, sev_obs = full_data_qlg.hs_class,
                sev_baseline = ynew_baseline_mode,
                qlg_treatment_type = full_data_qlg.qlg_treatment_type,
                prob_high_pred = prob_high, prob_high_baseline = prob_high_baseline)
CSV.write("../../data/qlg/pred.csv", out)

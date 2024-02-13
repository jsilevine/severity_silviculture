using CSV, DataFrames, MLJ, MLJDecisionTreeInterface,
    DecisionTree, Plots, StatsBase, Random, MLJBase,
    JLD2, Shapley, Shapley.ComputationalResources,
    NaNStatistics, Loess

## load in data
data = CSV.read("../data/complete_data.csv", DataFrame);
data = data[completecases(data),:];


nrow(data[data.own_type .== "Private Industrial" .&& data.hdw .< maximum(data[data.own_type .== "Private Industrial", :hdw]) .&& data.hs_class .== "high", :]) ./
    nrow(data[data.own_type .== "Private Industrial" .&& data.hdw .< maximum(data[data.own_type .== "Private Industrial", :hdw]), :])

nrow(data[data.own_type .== "Federal" .&& data.hdw .< maximum(data[data.own_type .== "Private Industrial", :hdw]) .&& data.hs_class .== "high", :]) ./
    nrow(data[data.own_type .== "Federal" .&& data.hdw .< maximum(data[data.own_type .== "Private Industrial", :hdw]), :])

mean(data[data.own_type .== "Private Industrial" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]) .&& data.hs_class .== "high", :hdw])
mean(data[data.own_type .== "Federal" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]) .&& data.hs_class .== "high", :hdw])




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

folds_formatted = load_object("../data/folds_formatted.jld2");

data.holdout .= 1
data[vcat(folds_formatted[1][1], folds_formatted[1][2]), :holdout] .= 0

data.hs_class .= "low-none";
data[data.CBI_bc .> 2.25, :hs_class] .= "high";
data[data.CBI_bc .< 2.25 .&& data.CBI_bc .> 1.25, :hs_class] .= "med";

## target
target = Vector(data[:,:hs_class]);
target = coerce(target, OrderedFactor);
levels!(target, ["low-none", "med", "high"]);

cols = [:clust, :mean_dens, :mean_area, :sd_area, :em,
        :elevation, :slope, :tpi, :heat_load, :prev_sev, :avg_fuel_moisture,
        :hdw, :avg_fuel_temp];
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
    CSV.write("../data/rf_param_search.csv", performance_data)
end

## check results
performance_data


##---------------------------------------------------------------
## Fit model with tuned parameters
##---------------------------------------------------------------

forest_model = MLJDecisionTreeInterface.RandomForestClassifier(n_trees = 40,
                                                               max_depth = 30,
                                                               n_subfeatures = 4)
fm = machine(forest_model,
             features,
             target)

## fit model, takes awhile
MLJ.fit!(fm, rows = vcat(folds_formatted[1][1], folds_formatted[1][2]))

save_object("../data/forest_model.jld2", fm)

##---------------------------------------------------------------
## check out of sample prediction performance and compute feature importances
##---------------------------------------------------------------
fm = load_object("../data/forest_model.jld2")

ho_list = findall(data.holdout .== 1)

yhat = MLJ.predict(fm, rows = ho_list)
yhat_mode = MLJ.predict_mode(fm, rows = ho_list)
yhat_mode = categorical(yhat_mode)
m(mode.(yhat), data[ho_list, :hs_class])
mcc(mode.(yhat), data[ho_list, :hs_class])

auc(yhat, data[ho_list, :hs_class])

mean(log_loss(yhat, data[ho_list, :hs_class]))

confusion_matrix(mode.(yhat), data[ho_list, :hs_class])

(92236 + 135317 + 443130) / length(yhat)


yprob = pdf(yhat, ["low-none", "med", "high"])
Plots.plot(data[ho_list, :mean_dens], yprob[:,3])

UnivariateFinite(['x', 'z'], [0.11111, 0.9], pool=missing, ordered=true)
##---------------------------------------------------------------
## Compute shapley metrics
##---------------------------------------------------------------
using Shapley: MonteCarlo

## create subset on which to compute shapley values
Random.seed!(2)
exp_list = sample(ho_list, 10000)
explain = copy(features[exp_list, :])

## compute shapley values
shap = shapley(η -> MLJ.predict(fm, η), MonteCarlo(CPUThreads(), 64), explain)

save_object("../data/shapley.jld2", shap)

## summarize shapley values
shapley_summary = DataFrame(feature = names(features),
                            mean_low_none = Vector{Float64}(undef, ncol(features)),
                            sd_low_none = Vector{Float64}(undef, ncol(features)),
                            mean_med = Vector{Float64}(undef, ncol(features)),
                            sd_med = Vector{Float64}(undef, ncol(features)),
                            mean_high = Vector{Float64}(undef, ncol(features)),
                            sd_high = Vector{Float64}(undef, ncol(features)))

for i in 1:length(shap)
    shp_mat = pdf(shap[i], ["low-none", "med", "high"])
    shapley_summary[i, [2,4,6]] .= vec(mean(abs.(shp_mat), dims = 1))
    shapley_summary[i, [3,5,7]] .= vec(std(abs.(shp_mat), dims = 1))
end

sort(shapley_summary, order(:mean_high, rev = true))

data_plot = sort(shapley_summary, order(:mean_high, rev = true))
bar(data_plot.mean_high, orientation = :h, yticks = (1:14, data_plot.feature),
    xlims = [0, 0.06], frame = :box, legend = :none, color = :black, title = "Feature importance (high)",
    yflip = true)
savefig("../plots/shapley_importance_high.pdf")

data_plot = sort(shapley_summary, order(:mean_med, rev = true))
bar(data_plot.mean_med, orientation = :h, yticks = (1:14, data_plot.feature),
    xlims = [0, 0.03], frame = :box, legend = :none, color = :black, title = "Feature importance (moderate)",
    yflip = true)
savefig("../plots/shapley_importance_med.pdf")

data_plot = sort(shapley_summary, order(:mean_low_none, rev = true))
bar(data_plot.mean_low_none, orientation = :h, yticks = (1:14, data_plot.feature),
    xlims = [0, 0.04], frame = :box, legend = :none, color = :black, title = "Feature importance (low-none)",
    yflip = true)
savefig("../plots/shapley_importance_low.pdf")

## summarize shapley values
shapley_summary = DataFrame(feature = names(features),
                            mean_low_none = Vector{Float64}(undef, ncol(features)),
                            sd_low_none = Vector{Float64}(undef, ncol(features)),
                            mean_med = Vector{Float64}(undef, ncol(features)),
                            sd_med = Vector{Float64}(undef, ncol(features)),
                            mean_high = Vector{Float64}(undef, ncol(features)),
                            sd_high = Vector{Float64}(undef, ncol(features)))

pi_sub = findall(data[exp_list, :own_type] .== "Private Industrial")
pf_sub = findall(data[exp_list, :own_type] .== "Federal")

length(shap)

for i in 1:length(shap)
    shp_mat = pdf(shap[i], ["low-none", "med", "high"])
    shapley_summary[i, [2,4,6]] .= vec(mean(abs.(shp_mat[pi_sub]), dims = 1))
    shapley_summary[i, [3,5,7]] .= vec(std(abs.(shp_mat[pi_sub]), dims = 1))
end

sort(shapley_summary, order(:mean_high, rev = true))

data_plot = sort(shapley_summary, order(:mean_high, rev = true))
bar(data_plot.mean_high, orientation = :h, yticks = (1:14, data_plot.feature),
    xlims = [0, 0.06], frame = :box, legend = :none, color = :black, title = "Feature importance (high -- Private)",
    yflip = true)
savefig("../plots/shapley_importance_high_private.pdf")


for i in 1:length(shap)
    shp_mat = pdf(shap[i], ["low-none", "med", "high"])
    shapley_summary[i, [2,4,6]] .= vec(mean(abs.(shp_mat[pf_sub]), dims = 1))
    shapley_summary[i, [3,5,7]] .= vec(std(abs.(shp_mat[pf_sub]), dims = 1))
end

sort(shapley_summary, order(:mean_high, rev = true))

data_plot = sort(shapley_summary, order(:mean_high, rev = true))
bar(data_plot.mean_high, orientation = :h, yticks = (1:14, data_plot.feature),
    xlims = [0, 0.06], frame = :box, legend = :none, color = :black, title = "Feature importance (high -- Public)",
    yflip = true)
savefig("../plots/shapley_importance_high_public.pdf")

Random.seed!(2)
pf_ss = sample(pf_sub, length(pi_sub))

balanced_sample = vcat(pf_ss, pi_sub)

for i in 1:length(shap)
    shp_mat = pdf(shap[i], ["low-none", "med", "high"])
    shapley_summary[i, [2,4,6]] .= vec(mean(abs.(shp_mat[pf_ss]), dims = 1))
    shapley_summary[i, [3,5,7]] .= vec(std(abs.(shp_mat[pf_ss]), dims = 1))
end

sort(shapley_summary, order(:mean_high, rev = true))

data_plot = sort(shapley_summary, order(:mean_high, rev = true))
bar(data_plot.mean_high, orientation = :h, yticks = (1:14, data_plot.feature),
    xlims = [0, 0.06], frame = :box, legend = :none, color = :black, title = "Feature importance (high -- Balanced)",
    yflip = true)
savefig("../plots/shapley_importance_high_balanced.pdf")

shp_mat = pdf(shap.mean_dens[pi_sub], ["low-none", "med", "high"])
Plots.plot(explain.mean_dens[pi_sub], shp_mat[:,3], seriestype = :scatter)

shp_mat = pdf(shap.clust, ["low-none", "med", "high"])
Plots.plot(explain.clust, shp_mat[:,3], seriestype = :scatter)

shp_mat = pdf(shap.hdw, ["low-none", "med", "high"])
Plots.plot(explain.hdw, shp_mat[:,3], seriestype = :scatter)

shp_mat = pdf(shap.avg_fuel_moisture, ["low-none", "med", "high"])
Plots.plot(explain.avg_fuel_moisture, shp_mat[:,3], seriestype = :scatter)

shp_mat = pdf(shap.em, ["low-none", "med", "high"])
Plots.plot(explain.em, shp_mat[:,3], seriestype = :scatter)

shp_mat = pdf(shap.mean_area, ["low-none", "med", "high"])
Plots.plot(explain.mean_area, shp_mat[:,3], seriestype = :scatter)

shp_mat = pdf(shap.sd_area, ["low-none", "med", "high"])
Plots.plot(explain.sd_area, shp_mat[:,3], seriestype = :scatter)


sum(data[data.own_type .== "Private Industrial", :hs]) / nrow(data[data.own_type .== "Private Industrial", :])
sum(data[data.own_type .== "Federal", :hs]) / nrow(data[data.own_type .== "Federal", :])

function hist_dif(data, col)
    dx = maximum(data[:, col]) / 150
    binedges = 0:dx:maximum(data[:, col])
    sp = copy(data[data.own_type .== "Private Industrial", col])
    sf = copy(data[data.own_type .== "Federal", col])
    aw = histcounts(sp[sample(1:length(sp), 500000)], binedges)
    bw = histcounts(sf[sample(1:length(sf), 500000)], binedges)
    bar(binedges, aw - bw, label = string(col) * " (Private - Federal)", bar_width = dx,
        frame = :box, color = :black)
end


hist_dif(data, :CBI_bc)

histogram(data.CBI_bc)
vline!([1.25, 2.25], legend = :none, title = "CBI (private - public)")

hist_dif(data[data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :], :CBI_bc)
vline!([1.25, 2.25], legend = :none, title = "CBI (private - public)")

hist_dif(data[data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :], :hdw)

string(:Stem_density)

nrow(data[data.hs_class .== "high" .&& data.own_type .== "Federal" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :]) /
    nrow(data[data.own_type .== "Federal" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :])

nrow(data[data.hs_class .== "med" .&& data.own_type .== "Federal" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :]) /
    nrow(data[data.own_type .== "Federal" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :])

nrow(data[data.hs_class .== "low-none" .&& data.own_type .== "Federal" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :]) /
    nrow(data[data.own_type .== "Federal" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :])

nrow(data[data.hs_class .== "high" .&& data.own_type .== "Private Industrial" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :]) /
    nrow(data[data.own_type .== "Private Industrial" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :])

nrow(data[data.hs_class .== "med" .&& data.own_type .== "Private Industrial" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :]) /
    nrow(data[data.own_type .== "Private Industrial" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :])

nrow(data[data.hs_class .== "low-none" .&& data.own_type .== "Private Industrial" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :]) /
    nrow(data[data.own_type .== "Private Industrial" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :])

nrow(data[data.hs_class .== "high" .&& data.own_type .== "Other" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :]) /
    nrow(data[data.own_type .== "Other" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :])

nrow(data[data.hs_class .== "med" .&& data.own_type .== "Other" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :]) /
    nrow(data[data.own_type .== "Other" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :])

nrow(data[data.hs_class .== "low-none" .&& data.own_type .== "Other" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :]) /
    nrow(data[data.own_type .== "Other" .&& data.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :])



hist_dif(data, :mean_dens)
savefig("../plots/mean_dens_comp.pdf")
hist_dif(data, :em)
savefig("../plots/em_comp.pdf")
hist_dif(data, :sd_area)
savefig("../plots/sd_area_comp.pdf")
hist_dif(data, :mean_area)
savefig("../plots/mean_area_comp.pdf")
hist_dif(data, :clust)
savefig("../plots/clust.pdf")
hist_dif(data, :elevation)
hist_dif(data, :hdw)

##---------------------------------------------------------------
## Partial dependency plots
##---------------------------------------------------------------

function pdp(;var = [:mean_dens], npts = 100, sp = 0.5, st = -1, var2_levels = [0], pts = true)

    bdf = DataFrame(fire_name = ["sugar", "walker", "north_complex", "sheep"],
                    clust = zeros(4),
                    mean_dens = zeros(4),
                    mean_area = zeros(4),
                    sd_area = zeros(4),
                    em = zeros(4),
                    elevation = zeros(4),
                    slope = zeros(4),
                    tpi = zeros(4),
                    heat_load = zeros(4),
                    prev_sev = zeros(4),
                    avg_fuel_moisture = zeros(4),
                    hdw = zeros(4),
                    avg_fuel_temp = zeros(4));
    varlist = [:fire_name, :clust, :mean_dens, :mean_area, :sd_area, :em, :elevation, :slope,
               :tpi, :heat_load, :prev_sev, :avg_fuel_moisture, :hdw, :avg_fuel_temp]

    if st == -1
        st = round((maximum(data[:,var[1]]) - minimum(data[:,var[1]])) / npts, digits = 2)
    end

    if length(var) > 1
        pd = DataFrame(Base.product(collect(range(0, round(maximum(data.mean_dens)), length = 100)),
                            var2_levels));
        rename!(pd, var)
        for v in varlist[2:14][varlist[2:14] .∉ [var]]
            pd[:,v] .= mean(data[:,v])
        end
    else
        pd = DataFrame()
        pd[:,var[1]] = collect(range(round(minimum(data[:,var[1]])), round(maximum(data[:,var[1]])), length = npts))
        for v in varlist[2:14][varlist[2:14] .!= var[1]]
            pd[:,v] .= mean(data[:,v])
        end
    end

    ## coerce to correct scitype
    for i in 1:ncol(pd)
        pd[:,i]  = coerce(float.(identity.(pd[:,i])), Continuous)
    end

    pd.fire_name .= "dixie";
    pd = vcat(pd, bdf);
    pd.fire_name = coerce(pd.fire_name, OrderedFactor);
    pd = pd[:, [:fire_name, :clust, :mean_dens, :mean_area, :sd_area, :em, :elevation, :slope,
                :tpi, :heat_load, :prev_sev, :avg_fuel_moisture, :hdw, :avg_fuel_temp]];

    yt = MLJ.predict(fm, pd);

    prs = pdf(yt, ["low-none", "med", "high"]);

    layout = @layout [a
                      b{1.0w,0.8h}]

    if length(var) > 1

        colors = ["#ffeda0", "#feb24c", "#f03b20"]

        np = npts * length(var2_levels)

        phigh = Plots.plot(layout = layout, link = :both, size = (500, 500), margin = -10Plots.px)
        if pts
            plot!(pd[1:np, var[1]], prs[1:np,3], subplot = 2, group = pd[1:np, var[2]],
                  seriestype = :scatter, markercolor = :match, color_palette = colors)
        end
        dx = maximum(data[:, var[1]]) / 150
        binedges = 0:dx:maximum(data[:, var[1]])
        spr = copy(data[data.own_type .== "Private Industrial", var[1]])
        sf = copy(data[data.own_type .== "Federal", var[1]])
        aw = histcounts(spr[sample(1:length(spr), 500000)], binedges)
        bw = histcounts(sf[sample(1:length(sf), 500000)], binedges)
        for i in 1:length(unique(pd[:,var[2]]))-1
            lm = loess(pd[pd[!,var[2]] .== unique(pd[!,var[2]])[i], var[1]], prs[pd[!,var[2]] .== unique(pd[!,var[2]])[i],3], span = 0.5)
            nx = range(extrema(pd[pd[!,var[2]] .== unique(pd[!,var[2]])[i], var[1]])...; step = 5)
            ny = Loess.predict(lm, nx)
            plot!(nx, ny, subplot = 2, color = colors[i], linewidth = 3, legend = :none, frame = :box)
        end
        xlabel!(string(var[1]))
        bar!(binedges, aw, subplot = 1, bar_width = dx, alpha = 0.3,
             frame = :none, color = :blue, legend = :none, linewidth = 0)
        bar!(binedges, bw, subplot = 1, bar_width = dx, alpha = 0.3,
             frame = :none, color = :green, legend = :none, linewidth = 0, title = "High")

        pmed = Plots.plot(layout = layout, link = :both, size = (500, 500), margin = -10Plots.px)
        if pts
            plot!(pd[1:np, var[1]], prs[1:np,2], group = pd[1:np, var[2]], subplot = 2,
                  seriestype = :scatter, markercolor = :match, color_palette = colors)
        end
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
             frame = :none, color = :green, legend = :none, linewidth = 0, title = "Med")

        plow = Plots.plot(layout = layout, link = :both, size = (500, 500), margin = -10Plots.px)
        if pts
            plot!(pd[1:np, var[1]], prs[1:np,1], group = pd[1:np, var[2]], subplot = 2,
                  seriestype = :scatter, markercolor = :match, color_palette = colors)
        end
        for i in 1:length(unique(pd[:,var[2]]))-1
            lm = loess(pd[pd[!,var[2]] .== unique(pd[!,var[2]])[i], var[1]], prs[pd[!,var[2]] .== unique(pd[!,var[2]])[i],1], span = 0.5)
            nx = range(extrema(pd[pd[!,var[2]] .== unique(pd[!,var[2]])[i], var[1]])...; step = 5)
            ny = Loess.predict(lm, nx)
            plot!(nx, ny, subplot = 2, color = colors[i], linewidth = 3, legend = :none, frame = :box)
        end
        xlabel!(string(var[1]))
        bar!(binedges, aw, subplot = 1, bar_width = dx, alpha = 0.3,
             frame = :none, color = :blue, legend = :none, linewidth = 0)
        bar!(binedges, bw, subplot = 1, bar_width = dx, alpha = 0.3,
             frame = :none, color = :green, legend = :none, linewidth = 0, title = "Low")

    else

        dx = maximum(data[:, var[1]]) / 150
        binedges = 0:dx:maximum(data[:, var[1]])
        spr = copy(data[data.own_type .== "Private Industrial", var[1]])
        sf = copy(data[data.own_type .== "Federal", var[1]])
        aw = histcounts(spr[sample(1:length(spr), 500000)], binedges)
        bw = histcounts(sf[sample(1:length(sf), 500000)], binedges)

        phigh = Plots.plot(layout = layout, link = :both, size = (400, 400), margin = -10Plots.px)
        if pts
            plot!(pd[1:npts, var[1]], prs[1:npts,3], seriestype = :scatter, color = :black, subplot = 2)
        end
        lm = loess(pd[1:npts, var[1]], prs[1:npts,3], span = sp)
        nx = range(extrema(pd[1:npts, var[1]])...; step = st)
        ny = Loess.predict(lm, nx)
        plot!(nx, ny, color = :red, linewidth = 3, legend = :none, frame = :box,
              ylims = [round(minimum(prs), digits = 2), round(maximum(prs), digits = 2)], subplot = 2)
        xlabel!(string(var[1]))
        bar!(binedges, aw, subplot = 1, bar_width = dx, alpha = 0.3,
             frame = :none, color = :blue, legend = :none, linewidth = 0)
        bar!(binedges, bw, subplot = 1, bar_width = dx, alpha = 0.3,
             frame = :none, color = :green, legend = :none, linewidth = 0, title = "High")

        pmed = Plots.plot(layout = layout, link = :both, size = (400, 400), margin = -10Plots.px)
        if pts
            plot!(pd[1:npts, var[1]], prs[1:npts,2], seriestype = :scatter, color = :black, subplot = 2)
        end
        lm = loess(pd[1:npts, var[1]], prs[1:npts,2], span = sp)
        nx = range(extrema(pd[1:npts, var[1]])...; step = st)
        ny = Loess.predict(lm, nx)
        plot!(nx, ny, color = :blue, linewidth = 3, legend = :none, frame = :box,
              ylims = [round(minimum(prs), digits = 2), round(maximum(prs), digits = 2)], subplot = 2)
        xlabel!(string(var[1]))
        bar!(binedges, aw, subplot = 1, bar_width = dx, alpha = 0.3,
             frame = :none, color = :blue, legend = :none, linewidth = 0)
        bar!(binedges, bw, subplot = 1, bar_width = dx, alpha = 0.3,
             frame = :none, color = :green, legend = :none, linewidth = 0, title = "Med")

        plow = Plots.plot(layout = layout, link = :both, size = (400, 400), margin = -10Plots.px)
        if pts
            plot!(pd[1:npts, var[1]], prs[1:npts,1], seriestype = :scatter, color = :black, subplot = 2)
        end
        lm = loess(pd[1:npts, var[1]], prs[1:npts,1], span = sp)
        nx = range(extrema(pd[1:npts, var[1]])...; step = st)
        ny = Loess.predict(lm, nx)
        plot!(nx, ny, color = :green, linewidth = 3, legend = :none, frame = :box,
              ylims = [round(minimum(prs), digits = 2), round(maximum(prs), digits = 2)], subplot = 2)
        xlabel!(string(var[1]))
        bar!(binedges, aw, subplot = 1, bar_width = dx, alpha = 0.3,
             frame = :none, color = :blue, legend = :none, linewidth = 0)
        bar!(binedges, bw, subplot = 1, bar_width = dx, alpha = 0.3,
             frame = :none, color = :green, legend = :none, linewidth = 0, title = "Low")

    end

    return Plots.plot(phigh, pmed, plow)

end

## forest structure
pdp(var = [:mean_dens], pts = false)
savefig("../plots/pdp_mean_dens.pdf")
pdp(var = [:clust], pts = false)
savefig("../plots/pdp_clust.pdf")
pdp(var = [:mean_area], pts = false)
savefig("../plots/pdp_mean_area.pdf")
pdp(var = [:sd_area], pts = false)
savefig("../plots/pdp_sd_area.pdf")
pdp(var = [:em], sp = 0.3, pts = false)
savefig("../plots/pdp_em.pdf")

## weather
pdp(var = [:hdw], pts = false)
savefig("../plots/pdp_hdw.pdf")
pdp(var = [:avg_fuel_moisture], sp = 0.3, pts = false)
savefig("../plots/pdp_fuel_moisture.pdf")
pdp(var = [:avg_fuel_temp], sp = 0.3, pts = false)
savefig("../plots/pdp_avg_fuel_temp.pdf")
pdp(var = [:prev_sev], pts = false)
savefig("../plots/pdp_prev_sev.pdf")

## topography
pdp(var = [:elevation], pts = false)
savefig("../plots/pdp_elevation.pdf")
pdp(var = [:slope], pts = false)
savefig("../plots/pdp_slope.pdf")
pdp(var = [:tpi], pts = false)
savefig("../plots/pdp_tpi.pdf")
pdp(var = [:heat_load], pts = false)
savefig("../plots/pdp_heat_load.pdf")

pdp(var = [:mean_dens, :hdw], var2_levels = [5, 75, 150], pts = false)
savefig("../plots/pdp_dens_hdw.pdf")

pdp(var = [:mean_dens, :em], var2_levels = [10, 15, 20], pts = false)

pdp(var = [:mean_dens, :prev_sev], var2_levels = [0.5, 1.5, 2.5], pts = false)
savefig("../plots/pdp_dens_prevsev.pdf")

pdp(var = [:hdw, :mean_dens], var2_levels = [50, 150, 250], pts = false)
savefig("../plots/pdp_hdw_dens.pdf")


##---------------------------------------------------------------
## predictions when only varying forest structure
##---------------------------------------------------------------

fake_features = copy(features);
fake_features.elevation .= mean(features.elevation);
fake_features.slope .= mean(features.slope);
fake_features.tpi .= mean(features.tpi);
fake_features.heat_load .= mean(features.heat_load);
fake_features.prev_sev .= mean(features.prev_sev);
fake_features.avg_fuel_moisture .= mean(features.avg_fuel_moisture);
fake_features.hdw .= mean(features.hdw);
fake_features.avg_fuel_temp .= mean(features.avg_fuel_temp);

fake_predict = MLJ.predict(fm, fake_features[ho_list, :]);

probs = pdf(yhat, ["low-none", "med", "high"]);

pi_list = findall(data[ho_list, :own_type] .== "Private Industrial" .&& data[ho_list, :elevation] .< 1500)
pf_list = findall(data[ho_list, :own_type] .== "Federal" .&& data[ho_list, :elevation] .< 1500)

mean(probs[pi_list,3])
mean(probs[pf_list,3])


mean(data[data.own_type .== "Private Industrial" .&& data.elevation .< 1500, :hdw])
mean(data[data.own_type .== "Federal" .&& data.elevation .< 1500, :hdw])


##---------------------------------------------------------------
## New model with balanced elevation
##---------------------------------------------------------------

forest_model = MLJDecisionTreeInterface.RandomForestClassifier(n_trees = 40,
                                                               max_depth = 30,
                                                               n_subfeatures = 4)

features_bal = copy(features[features.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :])
target_bal = copy(target[features.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]), :])

## coerce to correct scitype
for i in 1:length(cols)
    features_bal[:,cols[i]]  = coerce(float.(identity.(features_bal[:,cols[i]])), Continuous)
end

target_bal = coerce(target_bal, OrderedFactor);
scitype(features_bal[:,1])

fm_bal = machine(forest_model,
                 features,
                 target)


inlist = findall(features.elevation .< maximum(data[data.own_type .== "Private Industrial", :elevation]))
fitlist = vcat(folds_formatted[1][1], folds_formatted[1][2]) .∈ [inlist]
fl = vcat(folds_formatted[1][1], folds_formatted[1][2])[Fitlist]

ho_list


features_bal

## fit model, takes awhile
MLJ.fit!(fm_bal, rows = fl)

save_object("../data/forest_model_balanced.jld2", fm_bal)

ho_list = findall(data.holdout .== 1)
ho_list_bal = ho_list[ho_list .∈ [inlist]]

## create subset on which to compute shapley values
Random.seed!(2)
exp_list_bal = sample(ho_list_bal, 10000)
explain = copy(features[exp_list_bal, :])

## compute shapley values
shap_bal = shapley(η -> MLJ.predict(fm_bal, η), MonteCarlo(CPUThreads(), 64), explain)

save_object("../data/shapley_bal.jld2", shap_bal)

## summarize shapley values
shapley_summary = DataFrame(feature = names(features),
                            mean_low_none = Vector{Float64}(undef, ncol(features)),
                            sd_low_none = Vector{Float64}(undef, ncol(features)),
                            mean_med = Vector{Float64}(undef, ncol(features)),
                            sd_med = Vector{Float64}(undef, ncol(features)),
                            mean_high = Vector{Float64}(undef, ncol(features)),
                            sd_high = Vector{Float64}(undef, ncol(features)))

for i in 1:length(shap_bal)
    shp_mat = pdf(shap_bal[i], ["low-none", "med", "high"])
    shapley_summary[i, [2,4,6]] .= vec(mean(abs.(shp_mat), dims = 1))
    shapley_summary[i, [3,5,7]] .= vec(std(abs.(shp_mat), dims = 1))
end

sort(shapley_summary, order(:mean_high, rev = true))

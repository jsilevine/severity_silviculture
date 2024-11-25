cd("../../")

using Pkg

Pkg.activate("./")

using CSV, RData, DataFrames, Statistics, StatsBase, ProgressBars
using ArchGDAL; const AG = ArchGDAL
using Rasters, RasterDataSources, Distances
using Dates, Plots
using NearestNeighbors

## Need to combine:
## 1. Topography
## 2. canopy_cover (ladder_fuels)
## 3. cbi (severity)
## 4. daily progression data
## 5. weather
## 6. forest structure (lidar)


fire_names = ["dixie", "northcomplex", "sugar", "sheep", "walker"]

##---------------------------------------------------------------
## 1. Topography and Canopy Cover
##---------------------------------------------------------------

topography = CSV.read("data/Topography/derived/topography.csv", DataFrame);
select!(topography, [:x, :y, :elevation, :slope, :tpi, :heat_load]);

canopy_cover = CSV.read("data/canopy_cover/processed/cc_complete.csv", DataFrame);
select!(canopy_cover, ["x", "y", "cc_2.8", "cc_8.16", "cc_16.32", "cc_32.n"]);
rename!(canopy_cover, "cc_2.8" => :cc2_8, "cc_8.16" => :cc8_16, "cc_16.32" => :cc16_32, "cc_32.n" => :cc32_n);

full_data = leftjoin(canopy_cover, topography, on = [:x, :y]);

topography = nothing;
canopy_cover = nothing;

function calc_em(i, data = full_data)
    cc = collect(data[i, [:cc2_8, :cc8_16, :cc16_32, :cc32_n]])
    cc = cc[cc .> 0.0]
    cc_p = cc ./ sum(cc)
    em = ((1 / sum(cc_p.^2)) / length(cc)) * mean(cc)
    return em
end

full_data.em = Vector{Float64}(undef, nrow(full_data));
for i in ProgressBar(1:nrow(full_data))
    full_data[i, :em] = calc_em(i, full_data)
end

##---------------------------------------------------------------
## 2. Climate water deficit
##---------------------------------------------------------------
cwd = CSV.read("data/cwd/cwd.csv", DataFrame, select = [:x, :y, :cwd])

full_data = leftjoin(full_data, cwd, on = [:x, :y])
cwd = nothing;

##---------------------------------------------------------------
## 3. Severity
##---------------------------------------------------------------

for f in fire_names
    sev = CSV.read("data/cbi/" * f * "_processed.csv", DataFrame, select = [:x, :y, :CBI_bc])
    sev.fire_name .= f
    rename!(sev, :CBI_bc => :cbi)
    if f == fire_names[1]
        global severity = sev
    else
        global severity = vcat(severity, sev)
    end

end

full_data = leftjoin(severity, full_data, on = [:x, :y]);
full_data = full_data[.!ismissing.(full_data.cc2_8),:];

full_data.hs .= 0;
full_data[full_data.cbi .> 2.25, :hs] .= 1;
full_data.hs_class .= "not_high";
full_data[full_data.cbi .> 2.25, :hs_class] .= "high";

severity = nothing

##---------------------------------------------------------------
## 4. Daily progressions
##---------------------------------------------------------------

for f in fire_names
    db = CSV.read("data/fire_progression/smooth_round/" * f * "_progression.csv",
                  DataFrame, select = [:x, :y, :datetime])
    db.fire_name .= f
    db.datetime_formatted = round.(unix2datetime.(db.datetime .* 86400), Dates.Minute)

    if f == fire_names[1]
        global daily_burned = db
    else
        global daily_burned = vcat(daily_burned, db)
    end

end

full_data = leftjoin(full_data, daily_burned, on = [:x, :y, :fire_name]);
full_data = full_data[.!ismissing.(full_data.datetime_formatted), :];
daily_burned = nothing;


full_data.prev_sev = Vector{Float64}(undef, nrow(full_data));


function dwm_for_row(j, data1, data2, r = 10000)
    x_c = data2[j, :x]
    y_c = data2[j, :y]
    data1_sub = data1
    #data1_sub = data1[(data1.x .> x_c - r) .&& (data1.x .< x_c + r) .&&
    #                  (data1.y .> y_c - r) .&& (data1.y .< y_c + r), :]
    if nrow(data1_sub) > 0
        distdf = vcat(DataFrame(data2[j, [:x, :y]]), data1_sub[:, [:x, :y]])
        data1_sub.d .= pairwise(Euclidean(), eachrow(Array(distdf)))[2:nrow(data1_sub)+1,1]
        data1_sub = sort(data1_sub, :d)[1:100,:]
        weights = 1 ./ data1_sub.d
        weighted_sum = sum(weights .* data1_sub[:, :cbi])
        dwm_value = weighted_sum / sum(weights)
        return dwm_value
    else
        return 0.0
    end
end


function dwm_for_row(j, data1, data2, r = 10000)
    x_c = data2[j, :x]
    y_c = data2[j, :y]
    r = 10000

    data1_sub = data1
    #data1_sub = data1[(data1.x .> x_c - r) .&& (data1.x .< x_c + r) .&&
    #                  (data1.y .> y_c - r) .&& (data1.y .< y_c + r), :]
    #data1_sub = data1[sample(1:nrow(data1), 750), :]
    if nrow(data1_sub) == 0
        return 0.0
    else
        distdf = vcat(DataFrame(data2[j, [:x, :y]]), data1_sub[:, [:x, :y]]);
        d = pairwise(Euclidean(), eachrow(Array(distdf)))[:,1];
        weights = 1 ./ d[2:(size(data1_sub, 1) + 1), 1];
        weighted_sum = sum(weights .* data1_sub[:, :cbi]);
        dwm_value = weighted_sum / sum(weights);
        return dwm_value
    end
end


function dwm_for_date(i, fire, dt, full_data)
    data2 = full_data[full_data.fire_name .== fire .&& (full_data.datetime_formatted .== dt[i]), :]
    if i == 1
        return repeat([0.0], nrow(data2))
    else
        data1 = full_data[full_data.fire_name .== fire .&& (full_data.datetime_formatted .== dt[i - 1]), :]
        data1 = data1[sample(1:nrow(data1), 1000), :]
        output = Vector{Float64}(undef, nrow(data2))
        Threads.@threads for j in 1:nrow(data2)
            output[j] = dwm_for_row(j, data1, data2, 1000)
        end

        # Use `pmap` for parallel computation
        return output
    end
end


# Calculate weighted average severity on previous timestep
for fire in fire_names
    dt = unique(full_data[full_data.fire_name .== fire, :datetime_formatted])
    dt = sort(dt)
    for i in 1:length(dt)
        full_data[full_data.fire_name .== fire .&& (full_data.datetime_formatted .== dt[i]), :prev_sev] .= dwm_for_date(i, fire, dt, full_data)
        println("completed: ", fire, ", iteration ", i, " of ", length(dt))
    end
end


CSV.write("data/intermediate_test.csv", full_data)

##---------------------------------------------------------------
## 5. Weather
##---------------------------------------------------------------

cols_to_keep = [:x, :y, :datetime, :avg_air_temp, :max_air_temp, :avg_wind_speed, :max_wind_speed,
                :avg_relative_humidity, :avg_fuel_moisture, :avg_fuel_temp]

for f in fire_names
    w = CSV.read("data/weather/" * f * "_weather/complete_weather.csv",
                  DataFrame, select = cols_to_keep)
    if f == fire_names[1]
        global weather = w
    else
        global weather = vcat(weather, w)
    end
end

## convert R's datetime to something legible by Dates.jl
weather.datetime_formatted = Vector{DateTime}(undef, nrow(weather));
for i in ProgressBar(1:nrow(weather))
    date_str = weather[i, :datetime]
    cleaned_str = strip(date_str, ['(', ')'])
    weather[i, :datetime_formatted] = DateTime(cleaned_str, "mm/dd/yy HH:MM:SS") + Dates.Year(2000)
end

function calc_es(temp)
    6.11 * exp((2.5e6 / 461) * (1 / 273 - 1 / (273 + temp)))
end

function calc_vpd(rh, temp)
    ((100 - rh) / 100) * calc_es(temp)
end

weather.vpd = calc_vpd.(weather.avg_fuel_moisture, weather.avg_air_temp);
weather.hdw = weather.vpd .* weather.avg_wind_speed;

CSV.write("data/weather/full_weather.csv", weather);

weather_sub = select(weather, [:x, :y, :datetime_formatted, :hdw, :avg_fuel_moisture]);

full_data = leftjoin(full_data, weather_sub, on = [:x, :y, :datetime_formatted]);

weather = nothing;
weather_sub = nothing;

##---------------------------------------------------------------
## 6. Ownership
##---------------------------------------------------------------

ownership = CSV.read("data/ownership/ownership.csv", DataFrame,
                     select = [:x, :y, :own_int, :own_type])

full_data = leftjoin(full_data, ownership, on = [:x, :y])

##---------------------------------------------------------------
## 7. forest structure
##---------------------------------------------------------------

tree_data_30 = CSV.read("data/processed/tree_data_30m_full.csv", DataFrame,
                        select = [:x, :y, :clust, :mean_dens, :mean_ht, :median_ht, :max_ht]);

rename!(tree_data_30, :clust => :clust_30, :mean_dens => :mean_dens_30, :mean_ht => :mean_ht_30,
        :median_ht => :median_ht_30, :max_ht => :max_ht_30);

full_data = leftjoin(full_data, tree_data_30, on = [:x, :y]);
full_data = full_data[.!ismissing.(full_data.clust_30),:];

tree_data_180 = CSV.read("data/processed/tree_data_180m_full.csv", DataFrame,
                        select = [:x, :y, :clust, :mean_dens, :mean_ht, :median_ht, :max_ht]);

rename!(tree_data_180, :clust => :clust_180, :mean_dens => :mean_dens_180, :mean_ht => :mean_ht_180,
        :median_ht => :median_ht_180, :max_ht => :max_ht_180);

full_data = leftjoin(full_data, tree_data_180, on = [:x, :y]);

gap_data_30 = CSV.read("data/processed/gap_data_30m_full.csv", DataFrame,
                        select = [:x, :y, :mean_area, :median_area, :sd_area, :percent_open, :mean_frac]);

rename!(gap_data_30, :mean_area => :mean_area_30, :median_area => :median_area_30, :sd_area => :sd_area_30,
        :percent_open => :percent_open_30, :mean_frac => :mean_frac_30);

full_data = leftjoin(full_data, gap_data_30, on = [:x, :y]);
full_data = full_data[.!ismissing.(full_data.mean_area_30),:];

gap_data_180 = CSV.read("data/processed/gap_data_180m_full.csv", DataFrame,
                        select = [:x, :y, :mean_area, :median_area, :sd_area, :percent_open, :mean_frac]);

rename!(gap_data_180, :mean_area => :mean_area_180, :median_area => :median_area_180, :sd_area => :sd_area_180,
        :percent_open => :percent_open_180, :mean_frac => :mean_frac_180);

full_data = leftjoin(full_data, gap_data_180, on = [:x, :y]);
full_data = full_data[.!ismissing.(full_data.mean_area_180),:];

CSV.write("data/complete_data.csv", full_data)





##---------------------------------------------------------------
## SCRATCH
##---------------------------------------------------------------


using LinearAlgebra, Statistics, Plots

# Example: Point A and set of N points B
A = [1.0 2.0]  # Point A (x, y)
B = [3.0 2.0; 5.0 9.0; 4.0 8.0; 10.0 2.0]  # Set of N points B (each row is a point [x, y])

# Step 1: Compute the direction vectors from A to each point in B
vectors = B .- A  # Vectors from A to each point in B

# Step 2: Normalize each direction vector
norm_vectors = normalize.(eachrow(vectors))  # Normalize each direction vector

# Step 3: Compute the average direction vector
avg_direction = mean(norm_vectors, dims=1)  # Average the normalized direction vectors

# Step 4: Compute the slope (m_avg) of the average direction
x_avg, y_avg = avg_direction[1][1], avg_direction[1][2]
m_avg = y_avg / x_avg  # Slope of the average direction
b_avg = A[2] - m_avg * A[1]

# Result: The slope of the average direction
println("Average direction slope: ", m_avg)

m_perp = -1 / m_avg  # Slope of the perpendicular line
b_perp = A[2] - m_perp * A[1]  # Intercept of the perpendicular line

theta_perp = atan(m_perp)

# Step 3: Create the rotation matrix
rotation_matrix = [cos(theta_perp) -sin(theta_perp); sin(theta_perp) cos(theta_perp)]

# Step 4: Translate points in B so that A becomes the origin (subtract A from each point in B)
translated_points = B .- A  # Translate B points by subtracting A

# Step 5: Apply the rotation to each point in the translated points
rotated_points = (rotation_matrix * translated_points')'

# Step 6: Translate the points back by adding A
rotated_points = rotated_points .+ A  # Translate back by adding A

# Output the rotated points
println("Rotated points: ")
println(rotated_points)


scatter(vcat(A[1],B[:, 1]), vcat(A[2], B[:, 2]), label="Points B", legend=:none, ylim = [-10.0, 10.0], xlim = [-10.0, 10.0])
Plots.abline!(m_avg, b_avg)
Plots.abline!(m_perp, b_perp)
scatter!(rotated_points[:,1], rotated_points[:,2], label = "rotated")

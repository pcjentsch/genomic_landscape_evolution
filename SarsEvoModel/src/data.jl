using CSV, DataFrames
using Dates
using JSON
using DataStructures



datapath(fname) = joinpath(@__DIR__, "../data", fname)

using Serialization
function serial_load(load_func, fname)
    if !isfile(fname)
        @info "loading $fname"

        data = load_func()
        serialize(fname, data)
        return data
    else
        @info "loading serialized $fname..."

        return deserialize(fname)

    end
    @info "done"
end
function no_occurences(df)

    output_df = groupby(df, :lineage) |>
                df -> combine(nrow, df)

    output_df.pop_fraction = output_df.nrow ./ sum(output_df.nrow)
    sort!(output_df, :nrow)
    return output_df
end

struct Timeseries{T}
    dates::Vector{Date}
    values::Vector{T}
    function Timeseries(dates, values::AbstractVector{T}) where {T}
        (length(values) == length(dates)) || throw("arguments must be the same length!")
        return new{T}(dates, values)
    end
end
using Interpolations
function interpolate(xs::AbstractVector{Date}, ts::Timeseries{T}) where {T<:Number}
    begin_date = minimum(xs)
    @assert begin_date <= minimum(ts.dates)
    ts_xs = map(x -> Dates.value(x - begin_date), ts.dates)
    interp = LinearInterpolation(ts_xs, ts.values, extrapolation_bc=0.0)
    xs_vals = map(x -> Dates.value(x - begin_date), xs)
    return [interp(x) for x in xs_vals]
end

function interpolate(xs::AbstractVector{Date}, ts::Timeseries{T}) where {T<:AbstractArray{<:Real}}
    begin_date = minimum(xs)
    @assert begin_date <= minimum(ts.dates)
    ts_xs = map(x -> Dates.value(x - begin_date), ts.dates)
    interpolated = [zeros(size(ts.values[1])) for _ in xs]

    for i in 1:size(ts.values[1], 1), j in 1:size(ts.values[1], 2)
        vals_ij = [m[i, j] for m in ts.values]
        interp = LinearInterpolation(ts_xs, vals_ij, extrapolation_bc=0.0)
        xs_vals = map(x -> Dates.value(x - begin_date), xs)
        for (k, x) in enumerate(xs_vals)
            interpolated[k][i, j] = interp(x)
        end
    end
    return interpolated
end

struct LocationData
    dates::Vector{Date}
    stringency::Vector{Float64}
    cases_by_lineage::Vector{Matrix{Float64}}
    vaccination_mrna::Vector{Float64}
    function LocationData(stringency, vaccination, cases_by_lineage)

        first_date = min(
            minimum(stringency.dates),
            minimum(vaccination.dates),
            minimum(cases_by_lineage.dates),
        )

        last_date = max(
            maximum(stringency.dates),
            maximum(vaccination.dates),
            maximum(cases_by_lineage.dates),
        )
        dates = collect(first_date:Day(1):last_date)

        vaccination_interpolated = interpolate(dates, vaccination)
        stringency_interpolated = interpolate(dates, stringency)
        cases_by_lineage_interpolated = interpolate(dates, cases_by_lineage)


        return new(dates,
            stringency_interpolated,
            cases_by_lineage_interpolated,
            vaccination_interpolated)
    end
end

function OntarioLocationData()
    stringency = clean_strigency_data()
    vaccination = clean_vaccination_data()
    cases_by_lineage = load_cases_data()
    return LocationData(stringency, vaccination, cases_by_lineage)
end

function UKLocationData()
    country = "United Kingdom"
    owid_df = CSV.File(datapath("owid-covid-data.csv")) |> DataFrame |>
              df -> filter(:location => ==(country), df)

    vaccination_df = select(owid_df, [:date, :people_fully_vaccinated]) |>
                     df -> dropmissing(df)


    stringency_df = select(owid_df, [:date, :stringency_index]) |>
                    df -> dropmissing(df)

    cases_by_lineage = deserialize(datapath("uk_grids.data"))
    cases_by_lineage_ts = Timeseries(cases_by_lineage_dates, map(cases_by_lineage_to_matrix, cases_by_lineage))
    stringency_ts = Timeseries(stringency_df.date, stringency_df.stringency_index ./ 100)
    vaccination_ts = Timeseries(vaccination_df.date, vaccination_df.people_fully_vaccinated)

    return LocationData(stringency_ts, vaccination_ts, cases_by_lineage_ts)
end


function USALocationData()
    country = "United States"
    owid_df = CSV.File(datapath("owid-covid-data.csv")) |> DataFrame |>
              df -> filter(:location => ==(country), df)

    vaccination_df = select(owid_df, [:date, :people_fully_vaccinated]) |>
                     df -> dropmissing(df)


    stringency_df = select(owid_df, [:date, :stringency_index]) |>
                    df -> dropmissing(df)

    cases_df = select(owid_df, [:date, :new_cases_smoothed]) |>
               df -> dropmissing(df) |> df -> sort(df, :date)

    lineage_fractions, dates = deserialize(datapath("usa_grids.data"))
    cases_by_lineage = Matrix{Float64}[]
    cases_by_lineage_dates = Date[]
    for (lineage_fraction, date) in zip(lineage_fractions, dates)
        ind = findfirst(==(date), cases_df.date)
        if !isnothing(ind)
            push!(cases_by_lineage, lineage_fraction .* cases_df.new_cases_smoothed[ind])
            push!(cases_by_lineage_dates, date)
        end
    end
    cases_by_lineage_ts = Timeseries(cases_by_lineage_dates, cases_by_lineage)
    stringency_ts = Timeseries(stringency_df.date, stringency_df.stringency_index ./ 100)
    vaccination_ts = Timeseries(vaccination_df.date, vaccination_df.people_fully_vaccinated)

    return LocationData(stringency_ts, vaccination_ts, cases_by_lineage_ts)
end

function clean_vaccination_data()
    mrna_courses = [
        "Pfizer-BioNTech Comirnaty",
        "Pfizer-BioNTech Comirnaty pediatric 5-11 years",
        "mRNA mixed series",
        "Moderna Spikevax",
        "mRNA-protein subunit mixed series",
        "AstraZeneca Vaxzevria/COVISHIELD-mRNA mixed series",
        "Janssen mixed series",
    ]
    az_courses = [
        "Mixed series (unspecified)",
        "Janssen",
        "Novavax",
        "AstraZeneca Vaxzevria/COVISHIELD",
    ]
    vaccination_df = CSV.File(datapath("vaccination-coverage-byVaccineType.csv")) |> DataFrame
    vaccination_df.product_name = replace!(x -> x in mrna_courses ? "mrna" : "az", vaccination_df.product_name)
    cleaned_df = filter(:prename => ==("Ontario"), vaccination_df) |>
                 df -> groupby(df, [:week_end, :product_name]) |>
                       df -> combine(df, :numtotal_fully => sum)

    mrna_grouped = groupby(cleaned_df, :week_end) |> df -> combine(df, :numtotal_fully_sum => sum) |>
                                                           df -> sort(df, :week_end)
    mrna_timeseries = Timeseries(mrna_grouped.week_end, Float64.(mrna_grouped.numtotal_fully_sum_sum))
    return mrna_timeseries
end

function clean_strigency_data()
    stringency_df = CSV.File(datapath("COVID-19_STRINGENCY_INDEX.csv")) |> DataFrame
    return Timeseries(stringency_df.date, stringency_df.Ontario ./ 100.0)
end

function cases_by_lineage_to_matrix(cases_by_lineage::AbstractDataFrame)
    M = zeros(w, h)
    for (lineage, (x0, y0, width), date, pop, fraction) in eachrow(cases_by_lineage)
        for x in 1:w, y in 1:h
            M[y, x] += sigma(x - x0, y - y0; sigma_x=width, sigma_y=width, rounding=false) * pop
        end
    end
    return M
end

function load_cases_data()
    replace_dict = JSON.parsefile(datapath("lineages_replace.json"); dicttype=OrderedDict)
    lineage_data = CSV.File(datapath("on_lineage_report.csv")) |> DataFrame
    lineage_data.lineage = map(lineage_data.lineage) do lineage
        initial = first(split(lineage, "."))
        if haskey(replace_dict, initial)
            return replace(lineage, initial => replace_dict[initial])
        else
            return lineage
        end
    end

    covid_cases = CSV.File(datapath("covidtesting_on.csv")) |> DataFrame

    genomes_metadata = CSV.File(datapath("metadata_unzipped.tsv")) |> DataFrame
    genome_lineage = innerjoin(lineage_data, genomes_metadata; on=:taxon => :strain)
    genome_lineage.date = Date.(genome_lineage.date)
    select!(genome_lineage, [:lineage, :date])
    filter!(:lineage => !=("None"), genome_lineage)

    all_lineages = unique(genome_lineage.lineage)
    #sort by length, we want to match the longest tags preferentially
    tags = sort(first.(antigenic_map); by=length, rev=true)
    tag_per_lineage = map(all_lineages) do lineage
        ind = findfirst(occursin(lineage), tags)
        return isnothing(ind) ? nothing : tags[ind]
    end
    unmatched_lineages = [lineage for (tag, lineage) in zip(tag_per_lineage, all_lineages) if isnothing(tag)]
    @show unmatched_lineages
    lineage_tag_dict = Dict(tag => all_lineages[findall(==(tag), tag_per_lineage)] for tag in tags)


    lineage_fractions_by_date = DataFrame((date=first(df.date), lineages=no_occurences(df)) for df in groupby(genome_lineage, :date))
    cases_w_fractions = innerjoin(lineage_fractions_by_date, covid_cases; on=:date) |> df -> sort(df, :date)
    antigenic_landscape = map(antigenic_map) do (tag, (x, y, width))
        matching_lineages = lineage_tag_dict[tag]

        cases_by_date_for_tag = map(eachrow(cases_w_fractions)) do cases_on_date
            incident_cases = cases_on_date[Symbol("Confirmed Positive")]
            date = cases_on_date.date
            lineages_df = cases_on_date.lineages
            pop = 0.0
            fraction = 0
            for (; lineage, pop_fraction) in eachrow(lineages_df)
                if lineage in matching_lineages
                    pop += pop_fraction * incident_cases
                    fraction += pop_fraction
                end
            end
            return (; date, pop, fraction)
        end

        return (tag=tag, coords=(map_coords_to_model_space(x, y)..., 5 * width), cases=cases_by_date_for_tag)
    end |> DataFrame

    antigenic_landscape_flat = flatten(antigenic_landscape, :cases) |>
                               df -> transform(df, :cases => AsTable) |>
                                     df -> select(df, Not([:cases])) |>
                                           df -> groupby(df, :date)

    cases_by_lineage_dates = map(x -> x.date, keys(antigenic_landscape_flat))
    @assert issorted(cases_by_lineage_dates)
    cases_by_lineage = [antigenic_landscape_flat[k] for k in keys(cases_by_lineage_dates)]

    return Timeseries(cases_by_lineage_dates, map(cases_by_lineage_to_matrix, cases_by_lineage))
end
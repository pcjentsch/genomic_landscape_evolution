using CSV, DataFrames
using Dates
using JSON
using DataStructures



datapath(fname) = joinpath(@__DIR__, "../data", fname)

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
end

struct LocationData
    stringency::Timeseries{Float64}
    cases_by_lineage::Timeseries{Matrix{Float64}}
    vaccination_mrna::Timeseries{Float64}
    vaccination_az::Timeseries{Float64}
end

function load_covid_data(begin_date, end_date)

    replace_dict = JSON.parsefile(datapath("lineages_replace.json"); dicttype=OrderedDict)
    lineage_data = CSV.File(datapath("lineage_report.csv")) |> DataFrame
    lineage_data.lineage = map(lineage_data.lineage) do lineage
        initial = first(split(lineage, "."))
        if haskey(replace_dict, initial)
            return replace(lineage, initial => replace_dict[initial])
        else
            return lineage
        end
    end
    # display(lineage_data.lineage)

    covid_cases = CSV.File(datapath("covidtesting_on.csv")) |> DataFrame
    genomes_metadata = CSV.File(datapath("metadata_unzipped.tsv")) |> DataFrame
    genome_lineage = innerjoin(lineage_data, genomes_metadata; on=:taxon => :strain)
    select!(genome_lineage, [:lineage, :date])


    cases = filter(:date => x -> begin_date <= x <= end_date, covid_cases)
    lineage_fractions_by_date = DataFrame((date=Date(first(df).date), lineages=no_occurences(df)) for df in groupby(genome_lineage, :date))
    cases_w_fractions = innerjoin(lineage_fractions_by_date, cases; on=:date) |> df -> sort(df, :date)
    # before_begin_date = covid_cases[findlast(<=(begin_date), covid_cases[!, Symbol("Reported Date")]), Symbol("Total Cases")]
    # lineage_at_begin.populations = before_begin_date .* lineage_at_begin.pop_fraction


    antigenic_landscape = map(antigenic_map) do (tag, (x, y, width))
        cases_by_date_for_tag = map(eachrow(cases_w_fractions)) do cases_on_date
            incident_cases = cases_on_date[Symbol("Confirmed Positive")]
            date = cases_on_date.date
            lineages = cases_on_date.lineages
            pop = 0.0
            fraction = 0
            for (; lineage, pop_fraction) in eachrow(lineages)

                if occursin(tag, lineage)
                    display((lineage, tag))
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


    compare_cases = innerjoin(cases, combine(antigenic_landscape_flat, :pop => sum); on=:date)


    return cases_w_fractions, antigenic_landscape_flat, compare_cases
end
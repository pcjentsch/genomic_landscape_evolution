using CSV, DataFrames
using Dates
datapath(fname) = joinpath(@__DIR__, "../data", fname)

function no_occurences(df)

    output_df = groupby(df, :lineage) |>
                df -> combine(nrow, df)

    output_df.pop_fraction = output_df.nrow ./ sum(output_df.nrow)
    sort!(output_df, :nrow)
    return output_df
end
function load_covid_data(begin_date, end_date)
    lineage_data = CSV.File(datapath("lineage_report.csv")) |> DataFrame
    covid_cases = CSV.File(datapath("covidtesting_on.csv")) |> DataFrame
    genomes_metadata = CSV.File(datapath("subsampled_metadata_gisaid.tsv")) |> DataFrame
    genome_lineage = innerjoin(lineage_data, genomes_metadata; on = :taxon => :strain)
    select!(genome_lineage, [:lineage, :date_submitted])


    cases = filter(Symbol("Reported Date") => x -> begin_date <= x <= end_date, covid_cases)
    # genome_lineage.date .= trunc.(genome_lineage.date_submitted, Month)
    lineage_fractions_by_date = DataFrame((date = first(df).date_submitted, lineages = no_occurences(df)) for df in groupby(genome_lineage, :date_submitted))

    return lineage_fractions_by_date
    # before_begin_date = covid_cases[findlast(<=(begin_date), covid_cases[!, Symbol("Reported Date")]), Symbol("Total Cases")]
    # lineage_at_begin.populations = before_begin_date .* lineage_at_begin.pop_fraction


    # antigenic_landscape = map(antigenic_map) do (tag, (x, y, width))
    #     map(eachrow(cases)) do cases_on_date
    #         incident_cases  
    #         pop_at_tag = 0.0
    #         for (; lineage, populations) in eachrow(lineage_at_begin)
    #             if occursin(lineage, tag)
    #                 pop_at_tag += populations
    #             end
    #         end
    #     return (tag, (trunc(Int, x * 5), trunc(Int, y * 5), 5 * width), pop_at_tag)
    # end
end
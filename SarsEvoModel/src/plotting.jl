using Plots

const plotting_settings =
    (
        fontfamily="serif-roman",
        framestyle=:box,
        dpi=300,
        color_palette=palette(:seaborn_pastel)
    )


plots_path(fname; filetype="png") = joinpath(@__DIR__, "../plots", fname * "." * filetype)

function plot_antigenic_map()
    p = plot()
    for (name, (x, y, diameter)) in antigenic_map_paper
        scatter!(p, [x], [y];
            label=name, markersize=25 * diameter,
            size=(400, 300),
            legend=:outerright,
            xlims=(1, 9),
            ylims=(1, 7),
            xticks=collect(1:9),
            ytick=collect(1:7),
            xlabel="Antigenic distance",
            ylabel="Antigenic distance",
            plotting_settings...
        )
    end
    savefig(p, plots_path("antigenic_map_paper"))
end

function orfplot(recurrent_df)
    orf3a_coords = [
        ("Topological domain, Extracellular", Interval(1, 126)),
        ("Transmembrane, Helical", Interval(127, 183)),
        ("Topological domain, Cytoplasmic", Interval(184, 201)),
        ("Transmembrane, Helical", Interval(202, 279)),
        ("Topological domain, Extracellular", Interval(280, 303)),
        ("Transmembrane, Helical", Interval(304, 378)),
        ("Topological domain, Cytoplasmic", Interval(379, 825)),
    ]

    orf3a_df = groupby(filter(:gene_name => ==("ORF3a"), recurrent_df), :ind) |> df -> combine(df, :occurrence => sum)
    orf3a_df.description = map(eachrow(orf3a_df)) do r
        inds = filter(i -> (r.ind - 25393) in last(i), orf3a_coords)
        isempty(inds) && return missing
        return inds |> only |> first
    end
    dropmissing!(orf3a_df)
    orf3a_plot_scatter = plot()
    total_by_coords = map(orf3a_coords) do (label, coords)
        region_df = filter(:snp => s -> (s.ind - 25393) in coords, recurrent_df)
        return sum(region_df.occurrence) / (last(coords) - first(coords))
    end
    for (key, gdf) in pairs(groupby(orf3a_df, :description))

        scatter!(orf3a_plot_scatter, gdf.ind, gdf.occurrence_sum; label=key.description,
            ylabel="total recurrences", markersize=2, markerstrokewidth=0.1, plotting_settings...)
    end
    bar_colors = permutedims(plotting_settings.color_palette[[1, 2, 3, 2, 1, 2, 3]])
    bar_plot = bar(permutedims(string.(last.(orf3a_coords) .+ 25393)), permutedims(total_by_coords);
        color=bar_colors, legend=false, ylabel="total recurrences", xlabel="nucleotide index", plotting_settings...)

    savefig(plot(orf3a_plot_scatter, bar_plot; plotting_settings..., layout=(2, 1), size=(800, 500)), plots_path("orf3a"))
end

function plot_data(data::LocationData)
    incident_cases = data.cases_by_lineage
    anim = Animation()
    tlist = collect(1:length(data.dates))
    max_val = maximum(sum.(incident_cases))
    max_stringency = 1
    max_vac = maximum(data.vaccination_mrna)
    max_t = maximum(tlist)

    for (i, M) in enumerate(incident_cases)
        heatmap1 = heatmap(M; xlabel="antigenic distance", ylabel="antigenic distance", title="Incident Cases", seriescolor=cgrad(:Blues))
        ts1 = plot(tlist[1:i], sum.(incident_cases[1:i]); xlabel="time (days)", ylabel="total pop.", label="daily cases", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_val))
        ts2 = plot(tlist[1:i], data.vaccination_mrna[1:i]; xlabel="time (days)", ylabel="total pop.", label="total vaccinations", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_vac))
        ts3 = plot(tlist[1:i], sum.(data.stringency[1:i]); xlabel="time (days)", label="stringency", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_stringency))
        p = plot(heatmap1, ts1, ts2, ts3; layout=(4, 1), size=(400, 800),
            plotting_settings...)
        frame(anim, p)
    end
    gif(anim, plots_path("data"; filetype="gif"))
end
@views function plot_solution(sol, params)
    (; initial_population, location_data, begin_date) = params
    date_ind = findfirst(>=(begin_date), location_data.dates)
    incident_cases = location_data.cases_by_lineage[date_ind:end]
    initial_pop = initial_population * w * h
    anim = Animation()
    tlist = Float64[]
    S_ts = [u_t[1:w, 1:h] for u_t in sol.u]
    I_ts = [u_t[1:w, (h+1):(h*2)] for u_t in sol.u]

    C_ts = [u_t[1:w, (h*4+1):(h*5)] for u_t in sol.u]
    incident_ts = diff(sum.(C_ts))
    R_ts = [u_t[1:w, (h*2+1):(h*3)] for u_t in sol.u]
    V_ts = [u_t[1:w, (h*3+1):(h*4)] for u_t in sol.u]
    max_S = maximum(sum.(S_ts) ./ initial_pop)
    max_I = maximum(sum.(I_ts) ./ initial_pop)
    max_R = maximum(sum.(R_ts) ./ initial_pop)
    max_V = maximum(sum.(V_ts) ./ initial_population)
    max_incident = maximum(incident_ts ./ initial_population)
    tlist = sol.t
    display(length(tlist))
    max_t = maximum(tlist)
    # initial_population
    # p = plot(incident_ts)
    # plot!(p, sum.(incident_cases))
    # display(p)
    @showprogress for i in 1:length(tlist)-1

        heatmapS = heatmap(S_ts[i][2:end-1, 2:end-1]; xlabel="antigenic distance", ylabel="antigenic distance", title="S", seriescolor=cgrad(:Blues))
        heatmapI = heatmap(I_ts[i][2:end-1, 2:end-1]; xlabel="antigenic distance", title="I", seriescolor=cgrad(:Blues))
        # heatmapR = heatmap(R_ts[i][2:end-1, 2:end-1]; xlabel="antigenic distance", title="R", seriescolor=cgrad(:Blues))
        heatmapV = heatmap(V_ts[i][2:end-1, 2:end-1]; xlabel="antigenic distance", title="V", seriescolor=cgrad(:Blues))

        tsS = plot(tlist[1:i], sum.(S_ts[1:i]) ./ initial_pop; xlabel="time", ylabel="pop. fraction", label="susceptible", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_S))
        tsI = plot(tlist[1:i], incident_ts[1:i] ./ initial_population; xlabel="time", label="incident cases ", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_incident))
        # tsR = plot(tlist[1:i], sum.(R_ts[1:i]) ./ initial_pop; xlabel="time", label="recovered", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_R))
        tsV = plot(tlist[1:i], sum.(V_ts[1:i]) ./ initial_population; xlabel="time", label="vaccinated", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_V))
        tsData = plot!(tsI, tlist[1:i], sum.(incident_cases[1:i]) ./ initial_population;
            xlabel="time (days)", ylabel="total pop.", label="incident cases (data)",
            seriescolor=:Red, xlims=(0.0, max_t), ylims=(0.0, max(max_incident)))
        p = plot(
            heatmapS,
            heatmapI,
            # heatmapR,
            heatmapV,
            tsS,
            tsI,
            # tsR,
            tsV; margin=5Plots.mm, layout=(2, 3), size=(1200, 400),#, legend=:outerright,
            plotting_settings...)
        frame(anim, p)
    end
    gif(anim, plots_path("model"; filetype="gif"))
end

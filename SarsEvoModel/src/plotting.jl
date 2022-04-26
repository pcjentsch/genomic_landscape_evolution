using Plots
const plot_options = (
    dpi=300,
    fontfamily="Computer Modern",
    framestyle=:box
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
            plot_options...
        )
    end
    savefig(p, plots_path("antigenic_map_paper"))
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
        # display(M)
        heatmap1 = heatmap(M; xlabel="antigenic distance", ylabel="antigenic distance", title="Incident Cases", seriescolor=cgrad(:Blues))

        ts1 = plot(tlist[1:i], sum.(incident_cases[1:i]); xlabel="time (days)", ylabel="total pop.", label="daily cases", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_val))

        ts2 = plot(tlist[1:i], data.vaccination_mrna[1:i]; xlabel="time (days)", ylabel="total pop.", label="total vaccinations", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_vac))

        ts3 = plot(tlist[1:i], sum.(data.stringency[1:i]); xlabel="time (days)", label="stringency", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_stringency))

        p = plot(heatmap1, ts1, ts2, ts3; layout=(4, 1), size=(600, 1000),
            plot_options...)
        frame(anim, p)
    end
    gif(anim, plots_path("data"; filetype="gif"))
end
@views function plot_solution(sol)

    anim = Animation()
    tlist = Float64[]
    S_ts = [u_t[1:w, 1:h] for u_t in sol.u]
    I_ts = [u_t[1:w, (h+1):(h*2)] for u_t in sol.u]
    R_ts = [u_t[1:w, (h*2+1):(h*3)] for u_t in sol.u]
    V_ts = [u_t[1:w, (h*3+1):(h*4)] for u_t in sol.u]

    max_S = maximum(sum.(S_ts) ./ initial_pop)
    max_I = maximum(sum.(I_ts) ./ initial_pop)
    max_R = maximum(sum.(R_ts) ./ initial_pop)
    tlist = sol.t
    max_t = maximum(tlist)
    for (i, (S, I, R, V)) in enumerate(zip(S_ts, I_ts, R_ts, V_ts))

        heatmapS = heatmap(S; xlabel="antigenic distance", ylabel="antigenic distance", title="S", seriescolor=cgrad(:Blues))
        heatmapI = heatmap(I; xlabel="antigenic distance", title="I", seriescolor=cgrad(:Blues))
        heatmapR = heatmap(R; xlabel="antigenic distance", title="R", seriescolor=cgrad(:Blues))
        heatmapV = heatmap(V; xlabel="antigenic distance", title="V", seriescolor=cgrad(:Blues))

        tsS = plot(tlist[1:i], sum.(S_ts[1:i]) ./ initial_pop; xlabel="time", ylabel="pop. fraction", label="susceptible", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_S))
        tsI = plot(tlist[1:i], sum.(I_ts[1:i]) ./ initial_pop; xlabel="time", label="infected", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_I))
        tsR = plot(tlist[1:i], sum.(R_ts[1:i]) ./ initial_pop; xlabel="time", label="recovered", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_R))
        tsV = plot(tlist[1:i], sum.(V_ts[1:i]) ./ initial_pop; xlabel="time", label="vaccinated", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_R))

        p = plot(
            heatmapS,
            heatmapI,
            heatmapR,
            heatmapV,
            tsS,
            tsI,
            tsR,
            tsV;
            layout=(2, 4), size=(1200, 400),
            plot_options...)
        frame(anim, p)
    end
    gif(anim, plots_path("model"; filetype="gif"))
end

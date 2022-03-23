using Plots
const plot_options = (
    dpi=300,
    fontfamily="Computer Modern",
    framestyle=:box
)
plots_path(fname; filetype="png") = joinpath(@__DIR__, "../plots", fname * "." * filetype)

function plot_antigenic_map()
    p = plot()
    for (name, (x, y, diameter)) in antigenic_map
        scatter!(p, [x], [y];
            label=name, markersize=25 * diameter,
            size=(800, 600),
            legend=:outerright,
            xlims=(1, 9),
            ylims=(1, 7),
            xticks=collect(1:9),
            ytick=collect(1:7),
            xlabel="Antigenic distance",
            ylabel="Antigenic distance",
            title="Sars-CoV-2 antigenic map from Wilks et al (2022)",
            plot_options...
        )
    end
    savefig(p, plots_path("antigenic_map"))
end



function plot_solution(sol)

    anim = Animation()
    tlist = Float64[]
    S_ts = [u_t[1:w, 1:h] for u_t in sol.u]
    I_ts = [u_t[1:w, (h+1):(h*2)] for u_t in sol.u]
    R_ts = [u_t[1:w, (h*2+1):(h*3)] for u_t in sol.u]

    max_S = maximum(sum.(S_ts))
    max_I = maximum(sum.(I_ts))
    max_R = maximum(sum.(R_ts))
    tlist = sol.t
    max_t = maximum(tlist)
    for (i, (S, I, R)) in enumerate(zip(S_ts, I_ts, R_ts))

        heatmap1 = heatmap(S; xlabel="antigenic distance", ylabel="antigenic distance", title="S", seriescolor=cgrad(:Blues))
        heatmap2 = heatmap(I; xlabel="antigenic distance", title="I", seriescolor=cgrad(:Blues))
        heatmap3 = heatmap(R; xlabel="antigenic distance", title="R", seriescolor=cgrad(:Blues))

        ts1 = plot(tlist[1:i], sum.(S_ts[1:i]); xlabel="time", ylabel="total pop.", label="susceptible", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_S))
        ts2 = plot(tlist[1:i], sum.(I_ts[1:i]); xlabel="time", ylabel="total pop.", label="infected", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_I))
        ts3 = plot(tlist[1:i], sum.(R_ts[1:i]); xlabel="time", ylabel="total pop.", label="recovered", seriescolor=:Blue, xlims=(0.0, max_t), ylims=(0.0, max_R))

        p = plot(heatmap1, heatmap2, heatmap3, ts1, ts2, ts3;
            layout=(2, 3), size=(1000, 400),
            plot_options...)
        frame(anim, p)
    end
    gif(anim, plots_path("model"; filetype="gif"))
end

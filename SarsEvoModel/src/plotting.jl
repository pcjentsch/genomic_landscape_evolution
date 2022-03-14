using Plots
const plot_options = (
    dpi = 300,
    fontfamily = "Computer Modern",
    framestyle = :box
)
plots_path(fname; filetype = "png") = joinpath(@__DIR__, "../plots", fname * "." * filetype)

function plot_antigenic_map()
    p = plot()
    for (name, (x, y, diameter)) in antigenic_map
        scatter!(p, [x], [y];
            label = name, markersize = 25 * diameter,
            size = (800, 600),
            legend = :outerright,
            xlims = (1, 9),
            ylims = (1, 7),
            xticks = collect(1:9),
            ytick = collect(1:7),
            xlabel = "Antigenic distance",
            ylabel = "Antigenic distance",
            title = "Sars-CoV-2 antigenic map from Wilks et al (2022)",
            plot_options...
        )
    end
    savefig(p, plots_path("antigenic_map"))
end

function plot_kernel()
    p = heatmap(-15.0:15.0, -15.0:15.0, [sigma(i, j) for i in -15.0:15.0, j in -15.0:15.0]; plot_options...)
    savefig(p, plots_path("sigma"))
end



function plot_solution(sol)

    anim = Animation()
    tlist = Float64[]
    S_ts = Float64[]
    I_ts = Float64[]
    R_ts = Float64[]
    max_t = maximum(sol.t)
    for (i, (u_t, t)) in enumerate(zip(sol.u, sol.t))
        S = @view u_t[1:50, 1:50]
        I = @view u_t[1:50, 51:100]
        R = @view u_t[1:50, 101:150]
        push!(S_ts, sum(S))
        push!(I_ts, sum(I))
        push!(R_ts, sum(R))
        push!(tlist, t)

        heatmap1 = heatmap(S; xlabel = "antigenic distance", ylabel = "antigenic distance", title = "S", seriescolor = cgrad(:Blues))
        heatmap2 = heatmap(I; xlabel = "antigenic distance", title = "I", seriescolor = cgrad(:Blues))
        heatmap3 = heatmap(R; xlabel = "antigenic distance", title = "R", seriescolor = cgrad(:Blues))

        ts1 = plot(tlist, S_ts; xlabel = "time", ylabel = "total pop.", label = "S", seriescolor = :Blue, xlims = (0.0, max_t))
        ts2 = plot(tlist, I_ts; xlabel = "time", ylabel = "total pop.", label = "I", seriescolor = :Blue, xlims = (0.0, max_t))
        ts3 = plot(tlist, R_ts; xlabel = "time", ylabel = "total pop.", label = "R", seriescolor = :Blue, xlims = (0.0, max_t))

        p = plot(heatmap1, heatmap2, heatmap3, ts1, ts2, ts3;
            layout = (2, 3), size = (1500, 800),
            plot_options...)
        frame(anim, p)
    end
    gif(anim, plots_path("model"; filetype = "gif"))
end

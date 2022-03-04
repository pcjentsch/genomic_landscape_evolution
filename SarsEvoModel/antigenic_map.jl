using Plots

antigenic_map = [
    "P.1" => (4.2, 5.1, 0.4),
    "B.1.1.7" => (2.7, 3.8, 0.4),
    "D614G" => (2.7, 3.2, 0.4),
    "B.1.429" => (3.4, 2.2, 0.4),
    "B.1.617.2" => (2.25, 1.75, 0.4),
    "C.37" => (4.75, 2.25, 0.4),
    "B.1.617.1" => (5.6, 2.1, 0.4),
    "B.1.1.529" => (7.8, 1.5, 0.4),
    "B.1.526+E484K" => (5.6, 3.4, 0.4),
    "B.1.621" => (5.8, 4.7, 0.4),
    "B.1.351" => (5.3, 5.6, 0.4),
    "B.1.617.2(AY.2)+K417N" => (2.5, 1.75, 0.2),
    "B.1.617.2(AY.1)+K417N" => (2.2, 1.75, 0.2),
    "B.1.617.2+K417N" => (2.3, 2.0, 0.2),
    "B.1.617.2(AY.3)+E484Q" => (5.5, 1.8, 0.2),
    "B.1.1.7+E484K" => (5.2, 4.1, 0.2),
]

function plot_map()
    p = plot()
    for (name, (x, y, diameter)) in antigenic_map
        scatter!(p, [x], [y];
            label = name, markersize = 25 * diameter, dpi = 300, legend = :outerright,
            xlims = (0, 9),
            ylims = (0, 7),
            xticks = collect(0:9),
            ytick = collect(0:7),
            xlabel = "Antigenic distance",
            ylabel = "Antigenic distance")
    end
    savefig(p, "antigenic_map.png")
end
plot_map()
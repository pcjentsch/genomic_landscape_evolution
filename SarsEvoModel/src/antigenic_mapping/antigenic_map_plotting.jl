

function plot_mds(genomes_df, unique_df, lineages)
    # df = innerjoin(genomes_df, unique_df; on=:filtered_genome => :filtered_genome) |> df ->
    #     innerjoin(lineages, df; on=:id)
    # sort!(df, :date)
    unique_df = innerjoin(unique_df, lineages; on=:id => :taxon, makeunique=true)

    p = plot()
    display(unique_df)
    # unique_df.is_omicron = map(s -> occursin("BA", s), unique_df.Pangolin)
    for (key, gdf) in pairs(groupby(unique_df, :scorpio_call))
        p = scatter!(p, gdf.mds_x, gdf.mds_y; label=key.scorpio_call, xlabel="MDS1", ylabel="MDS2",
            title="Approximate antigenic cartography w/ RBD, UK samples",
            markersize=1.5,
            markerstrokewidth=0.3,
            size=(800, 500),
            plotting_settings...)
    end
    savefig(p, "multidimensional_scaling.png")

    # anim = Animation()
    # xlims = extrema(df.mds_x)
    # ylims = extrema(df.mds_y)
    # mds_x_accum = Float64[]
    # mds_y_accum = Float64[]
    # for (i, (key, gdf)) in enumerate(pairs(groupby(df, :date)))
    #     display(i)
    #     pts_df = unique(gdf, [:mds_x, :mds_y])
    #     append!(mds_x_accum, pts_df.mds_x)
    #     append!(mds_y_accum, pts_df.mds_y)
    #     p = scatter(mds_x_accum, mds_y_accum; label="samples", xlabel="MDS1", ylabel="MDS2",
    #         title="Estimated antigenic distance, $(key.date)",
    #         markersize=1.5,
    #         markerstrokewidth=0.3,
    #         xlims=xlims,
    #         ylims=ylims,
    #         plotting_settings...)
    #     frame(anim, p)
    # end
    # gif(anim, "animation.gif")
end


function orfplot()
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
    orf3a_plot = plot()
    for (key, gdf) in pairs(groupby(orf3a_df, :description))
        scatter!(orf3a_plot, gdf.ind, gdf.occurrence_sum; label=key.description,
            xlabel="nucleotide index", ylabel="total recurrences",
            title="orf3a recurrences by index", markersize=2, markerstrokewidth=0.1, plotting_settings...)
    end
    savefig(orf3a_plot, "orf3a_scatter_2.png")
end

function plot_samples(genomes_w_metadata)
    df_nonempty = filter(:in_rbd => g -> !isempty(g), genomes_w_metadata)
    datelist = [count(==(d), df_nonempty.date) / count(==(d), genomes_w_metadata.date) for d in sort(unique(df_nonempty.date))]
    numlist = [count(==(d), df_nonempty.date) for d in sort(unique(df_nonempty.date))]
    p = plot(sort(unique(df_nonempty.date)), datelist; legend=false, xlabel="date", ylabel="fraction of samples", title="Fraction of UK samples with mutations in RBD", plotting_settings...)
    savefig(p, "rbd_fraction.png")
    p = plot(sort(unique(df_nonempty.date)), numlist; legend=false, xlabel="date", ylabel="no samples", title="number of UK samples by date", plotting_settings...)
    savefig(p, "num_samples.png")
    df = groupby(genomes_w_metadata, :date) |>
         df -> combine(df,
        :binding_retained => mean => :mean,
        :binding_retained => std => :std,
        :binding_retained => (x -> maximum(x)) => :upper_quartile,
        :binding_retained => (x -> minimum(x)) => :lower_quartile
    ) |> df -> sort(df, :date)

    p = plot(df.date, df.mean; ribbon=df.std)
end

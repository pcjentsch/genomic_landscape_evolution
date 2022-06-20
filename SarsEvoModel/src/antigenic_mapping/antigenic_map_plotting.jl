

function plot_mds(unique_df)
    # df = innerjoin(genomes_df, unique_df; on=:filtered_genome => :filtered_genome) |> df ->
    #     innerjoin(lineages, df; on=:id)
    # sort!(df, :date)

    p = plot()
    display(unique_df)
    # unique_df.is_omicron = map(s -> occursin("BA", s), unique_df.Pangolin)
    for (key, gdf) in pairs(groupby(unique_df, :closest_mapped_lineage))
        p = scatter!(p, gdf.mds_x, gdf.mds_y; label=key.closest_mapped_lineage, xlabel="MDS1", ylabel="MDS2",
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

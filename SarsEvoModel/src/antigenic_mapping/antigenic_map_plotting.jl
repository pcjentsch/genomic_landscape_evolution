
#take 2 genomes, find common ancestor, find closest sampled genome to that ancestor

function plot_mds(fname, unique_df, dm)
    # df = innerjoin(genomes_df, unique_df; on=:filtered_genome => :filtered_genome) |> df ->
    #     innerjoin(lineages, df; on=:id)
    # sort!(df, :date)
    p = plot()
    for (key, gdf) in pairs(groupby(unique_df, :closest_mapped_lineage))
        p = scatter!(p, gdf.mds_x, gdf.mds_y; label=key.closest_mapped_lineage, xlabel="MDS1", ylabel="MDS2",
            markersize=1.5,
            markerstrokewidth=0.3,
            margin=5Plots.mm,
            size=(800, 300),
            legend=:outerright,
            plotting_settings...)
    end
    savefig(p, plots_path("$(fname)_multidimensional_scaling"))
    # mstress = 6
    # stress_list = zeros(mstress)
    # for i in 1:mstress
    #     mds = fit(MDS, dm; distances=true, maxoutdim=i)
    #     stress_list[i] = stress(mds)
    # end
    # p = plot(1:mstress, stress_list; xlabel="MDS Out-dimension", ylabel="Stress", plotting_settings...)
    # savefig(p, plots_path("$(fname)_mds_stress"))
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

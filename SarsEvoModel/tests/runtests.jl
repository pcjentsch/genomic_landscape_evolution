using SarsEvoModel

unique_df, filter_df, snp_weight_dict = SarsEvoModel.serial_load(
    false,#SarsEvoModel.unique_genomes(genomes_w_metadata, recurrent_df, binding_sites),
    SarsEvoModel.datapath("unique_df_usa.data")
)

using BenchmarkTools
@btime SarsEvoModel.dist(filter_df.filtered_genome[1], filter_df.filtered_genome[3], snp_weight_dict)
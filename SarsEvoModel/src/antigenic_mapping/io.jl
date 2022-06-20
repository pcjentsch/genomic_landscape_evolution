

function get_data_fasta(fasta_path, metadata_path, binding_sites)
    snps = SNPs_from_fastas(fasta_path) |> DataFrame
    metadata = CSV.File(metadata_path) |> DataFrame
    filter!(:Collection_Date => d -> d ∉ ("2020", "2021"), metadata)

    genomes_w_metadata = innerjoin(snps, metadata; on=:id => :strain)
    genomes_w_metadata.date = Date.(genomes_w_metadata.Collection_Date)
    genomes_w_metadata.in_rbd = map(genomes_w_metadata.genome) do snps
        return Set(snp for snp in snps if snp.ind in binding_sites)
    end
    unique_genomes = select(unique(genomes_w_metadata, :in_rbd), [:in_rbd])
    unique_genomes.binding_retained = @showprogress map(compute_binding_retained, unique_genomes.in_rbd)
    genomes_w_metadata = innerjoin(unique_genomes, genomes_w_metadata; on=:in_rbd)
    lineages_df = load_lineages()

    genomes_w_metadata = innerjoin(lineages_df, genomes_w_metadata; on=:taxon => :id)
    replace_dict = JSON.parsefile("SarsEvoModel//data/lineages_replace.json"; dicttype=OrderedDict)

    genomes_w_metadata.lineage = map(genomes_w_metadata.lineage) do lineage
        initial = first(split(lineage, "."))
        if haskey(replace_dict, initial)
            return replace(lineage, initial => replace_dict[initial])
        else
            return lineage
        end
    end
    omicron_prefix = "B.1.1.529"
    delta_prefix = "B.1.617.2"
    alpha_prefix = "B.1.1.7"

    genomes_w_metadata.category = map(genomes_w_metadata.lineage) do lineage
        occursin(omicron_prefix, lineage) && return 0 #"omicron"
        occursin(delta_prefix, lineage) && return 1 #"delta"
        occursin(alpha_prefix, lineage) && return 2 #"alpha"
        return 3 #"other"
    end
    filter!(∉(("Unassigned", "unclassifiable")), genomes_w_metadata)



    return genomes_w_metadata
end

function parse_recurrent_mutations(fpath)
    recurrent = CSV.File(fpath) |> DataFrame
    recurrent.snp = map(get_snp, recurrent.ID)
    recurrent.freq = recurrent.occurrence ./ sum(recurrent.occurrence)
    recurrent.ref = map(first, recurrent.ID)
    sort!(recurrent, :freq; rev=true)
    name_by_ind(ind) = gene_names[findall(interval -> ind in interval, gene_index)]
    recurrent.ind = map(snp -> snp.ind, recurrent.snp)
    recurrent.gene_name = map(name_by_ind, recurrent.ind)
    filter!(:gene_name => n -> length(n) == 1, recurrent)
    recurrent.gene_name = only.(recurrent.gene_name)
    return recurrent
end

function SNPs_from_record(sample_record, reference)
    undef_base = gap(DNA)
    id = FASTA.identifier(sample_record)
    genome = SNP[]
    sample_sequence = FASTA.sequence(sample_record)
    for (ind, (sample_base, ref_base)) in enumerate(zip(sample_sequence, reference))
        if ref_base != undef_base && sample_base != undef_base && ref_base != sample_base
            snp = SNP(ind, sample_base)
            push!(genome, snp)
        end
    end
    return (; id, genome)
end


function SNPs_from_fastas(aligned_fastas_dir::String)
    reference_reader = FASTA.Reader(open(joinpath(@__DIR__, "../../data/reference.fasta")))
    reference = FASTA.sequence(only(reference_reader))
    files = readdir(aligned_fastas_dir) |>
            l -> filter(s -> occursin("aligned", s), l) |>
                 l -> map(s -> joinpath(aligned_fastas_dir, s), l)
    display(files)
    file_streams = map(FASTA.Reader ∘ GzipDecompressorStream ∘ open, files)
    SNPs_by_sample = ThreadsX.map(r -> SNPs_from_record(r, reference), Iterators.flatten(file_streams))
    foreach(close, file_streams)
    return SNPs_by_sample
end



function load_lineages()
    lineages = CSV.File(datapath("uk_lineages.csv")) |> DataFrame
    lineages.id = replace.(lineages.taxon, "/" => "_", "-" => "N")
    lineages.scorpio_call_simplified = map(lineages.scorpio_call) do scorpio_call
        ismissing(scorpio_call) && return missing
        if occursin("Delta", scorpio_call)
            return "delta"
        else
            return "everything else"
        end
    end
    return lineages
end

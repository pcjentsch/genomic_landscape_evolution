

function get_data_fasta(fasta_path, metadata_path, binding_sites, lineage_path)
    if isfile("genomes_w_metadata.arrow")
        genomes_w_metadata = Arrow.Table("genomes_w_metadata.arrow") |> DataFrame
    else
        genomes_w_metadata = SNPs_from_fastas(fasta_path,metadata_path,binding_sites) |> DataFrame
        Arrow.write("genomes_w_metadata.arrow", genomes_w_metadata)
    end
    unique_genomes = select(unique(genomes_w_metadata, :all_in_rbd), [:all_in_rbd])
    unique_genomes.binding_retained = @showprogress map(compute_binding_retained, unique_genomes.all_in_rbd)
    genomes_w_metadata = innerjoin(genomes_w_metadata, unique_genomes; on = :all_in_rbd)
    display(genomes_w_metadata)
    lineages_df = load_lineages(lineage_path)
    genomes_w_metadata = innerjoin(lineages_df, genomes_w_metadata; on=:taxon => :all_ids)
    replace_dict = JSON.parsefile(joinpath(@__DIR__, "../../data/lineages_replace.json"); dicttype=OrderedDict)

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
    display(names(genomes_w_metadata))
    rename!(genomes_w_metadata, :all_in_rbd => :in_rbd)
    rename!(genomes_w_metadata, :all_genomes => :genome)
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
    undef_base = '-'
    id = FASTA.identifier(sample_record)
    genome = SNP[]
    sample_sequence = FASTA.sequence(sample_record)

    for (ind, (sample_base, ref_base)) in enumerate(zip(sample_sequence, reference))
        if ref_base != undef_base && sample_base != undef_base && ref_base != sample_base
            snp = SNP(ind, sample_base)
            push!(genome, snp)
        end
    end
    return (id, genome)
end
ArrowTypes.arrowname(::Type{SNP}) = :SNP
ArrowTypes.JuliaType(::Val{:SNP}) = SNP

function SNPs_from_fastas(aligned_fastas_dir::String, metadata_path,binding_sites)
    reference_reader = FASTA.Reader(open(joinpath(@__DIR__, "../../data/reference.fasta")))
    reference = FASTA.sequence(only(reference_reader))
    files = readdir(aligned_fastas_dir) |>
            l -> filter(s -> occursin("aligned", s), l) |>
                 l -> map(s -> joinpath(aligned_fastas_dir, s), l)

    metadata = CSV.File(metadata_path) |> DataFrame
    filter!(:Collection_Date => d -> d ∉ ("2020", "2021"), metadata)
    metadata.date = Date.(metadata.Collection_Date)
    all_genomes = Vector{Vector{SNP}}()
    all_ids = Vector{String}()
    all_in_rbd = Vector{Set{SNP}}()

    for (i,file) in enumerate(files)
        FASTAReader(GzipDecompressorStream(open(file))) do reader
            records = collect(reader)
            genomes = Vector{Vector{SNP}}(undef, length(records)) 
            ids = Vector{String}(undef, length(records))
            display("parsing file $file")
            Threads.@threads for j in 1:length(records)
                record = records[j]
                id, genome = SNPs_from_record(record, reference) 
                ids[j] = id
                genomes[j] = genome
            end
            in_rbd = map(g -> Set(snp for snp in g if snp.ind in binding_sites), genomes)
            display("done")

            append!(all_genomes, genomes)
            append!(all_ids, ids)
            append!(all_in_rbd, in_rbd)
        end
    end
    genomes_df = DataFrame((;all_ids, all_genomes, all_in_rbd))
    genomes_w_metadata = innerjoin(genomes_df, metadata; on=:all_ids => :name)
    return genomes_w_metadata
end



function load_lineages(path)
    lineages = CSV.File(path) |> DataFrame
    lineages.id = replace.(lineages.taxon, "/" => "_", "-" => "N")
    lineages.scorpio_call_simplified = map(lineages.scorpio_call) do scorpio_call
        ismissing(scorpio_call) && return missing
        if occursin("Delta", scorpio_call)
            return "delta"
        else
            return "everything else"
        end
    end
    if all(occursin.(" ", lineages.taxon))
        lineages.taxon = map(lineages.taxon) do id
            return split(id, " ") |> first
        end
    end
    return lineages
end

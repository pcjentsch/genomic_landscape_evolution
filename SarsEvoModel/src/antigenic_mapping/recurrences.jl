using BioSequences
function format_ref_w_snps(snp_list)
    reference_reader = FASTA.Reader(open(joinpath(@__DIR__, "../data/wuhCor1.fa")))
    reference = FASTA.sequence(only(reference_reader))[1:end-2]
    mutated_seq = deepcopy(reference)
    gene_starts = first.(sort(gene_index))
    for snp in snp_list
        (; ind, alt) = snp
        display(mutated_seq[ind])
        mutated_seq[ind] = alt
        display(mutated_seq[ind])
    end

    ref_aa = BioSequences.translate(reference)
    mutated_aa = BioSequences.translate(mutated_seq)# code=standard_genetic_code, allow_ambiguous_codons=true)
    for (; ind, alt) in snp_list
        gene_start = gene_starts[findlast(<=(ind), gene_starts)]
        aa_ind = floor(Int, ind / 3)
        display((reference[aa_ind*3:aa_ind*3+2], mutated_seq[aa_ind*3:aa_ind*3+2]))
        gene_aa_ind = floor(Int, (ind - gene_start) / 3)
        display((ref_aa[aa_ind], gene_aa_ind, mutated_aa[aa_ind]))
    end
    # do this biopython
    # grab coordinates
    # make vcf
    # seq_writer = FASTA.Writer(open(joinpath(@__DIR__, "../data/mutated_ref.fasta"), "w"))
    # rec = FASTA.Record("mutated", reference)
    # write(seq_writer, rec)
    # close(seq_writer)
    # close(reference_reader)
end

function analyze_recurrences()
    recurrent_df = parse_recurrent_mutations(joinpath(@__DIR__, "../data/filtered_mutations.tsv"))
    reference_reader = FASTA.Reader(open(joinpath(@__DIR__, "../data/reference.fasta")))
    reference = FASTA.sequence(only(reference_reader))[1:end-2]
    n = 0
    df1 = (filter([:snp, :ref] => (s, r) -> reference[s.ind] != convert(DNA, r), recurrent_df))
    recurrent_df.actual_ref = map(s -> reference[s.ind], recurrent_df.snp)
    df2 = filter([:ref, :actual_ref] => (r, ra) -> ra != convert(DNA, r), recurrent_df)
    return df1, df2
end

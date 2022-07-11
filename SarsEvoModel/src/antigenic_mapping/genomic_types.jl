struct SNP
    ind::UInt16
    alt::Char
end


function Base.isequal(a::SNP, b::SNP)
    return a.ind == b.ind && a.alt == b.alt
end

function Base.hash(a::SNP)
    return hash((a.ind, a.alt))
end


function weighted_dist(g1, g2, SNP_weights)
    d = 0
    for snp in g1
        if snp ∉ g2
            d += SNP_weights[snp]
        end
    end
    for snp in g2
        if snp ∉ g1
            d += SNP_weights[snp]
        end
    end
    return d
end

function snp_dist(g1, g2)
    d = 0
    for snp in g1
        if snp ∉ g2
            d += 1#SNP_weights[snp]
        end
    end
    for snp in g2
        if snp ∉ g1
            d += 1#SNP_weights[snp]
        end
    end
    return d
end





# function sorted_hash_dist(g1, g2)
#     d = 0
#     ind_1 = 1
#     ind_2 = 1
#     while ind
#     for e in g1
#         if e ∉ g2
#             d += 1
#         end
#     end
#     for e in g2
#         if e ∉ g1
#             d += 1
#         end
#     end
#     return d
# end
function get_snp(line::T) where {T<:AbstractString}
    alt = line[end]
    ind = parse(UInt16, line[2:end-1])
    return SNP(ind, alt)
end

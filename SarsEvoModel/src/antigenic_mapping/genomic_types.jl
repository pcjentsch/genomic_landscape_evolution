struct SNP
    ind::UInt16
    alt::Char
end


function Base.isequal(a::SNP, b::SNP)
    return a.ind == b.ind && a.ref == b.ref && a.alt == b.alt
end

function Base.hash(a::SNP)
    return hash((a.ind, a.ref, a.alt))
end


function dist(g1::Set{SNP}, g2::Set{SNP}, SNP_weights)
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
function get_snp(line::T) where {T<:AbstractString}
    alt = line[end]
    ind = parse(UInt16, line[2:end-1])
    return SNP(ind, alt)
end

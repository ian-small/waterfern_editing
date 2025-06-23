using CSV, DataFrames, StatsBase, BioSequences, FASTX, BioAlignments

const dna_stops = LongDNA{4}.([[DNA_T, DNA_A, DNA_A], [DNA_T, DNA_G, DNA_A], [DNA_T, DNA_A, DNA_G]])
const rna_stops = convert.(LongRNA{4}, dna_stops)

struct Feature
    contig::String
    type::String
    start::Int
    stop::Int
    strand::Char
    phase::Int
    attributes::Dict{String, String}
end

struct Gene
    id::String
    gene::Feature
    mRNAs::Vector{Feature}
    tRNAs::Vector{Feature}
    rRNAs::Vector{Feature}
    ncRNAs::Vector{Feature}
    CDSs::Vector{Feature}
    exons::Vector{Feature}
    introns::Vector{Feature}
    esites::Vector{Feature}
end

function attributes2dict(attribute_string::String)
    dict = Dict{String, String}()
    attributes = split(attribute_string, ";")
    for attribute in attributes
        kvs = split(attribute, "=")
        dict[first(kvs)] = last(kvs)
    end
    dict
end

function createFeature(row)
    Feature(row.sequence, row.feature, row.start, row.stop, row.strand[1], row.phase == "." ? 0 : parse(Int, row.phase), attributes2dict(row.attributes))
end

function add_feature(gene, feature)
    #check feature is on same contig as the gene
    @assert feature.contig == gene.gene.contig "contig error: $(string(feature)) $(string(gene.gene))"
    #check feature is on same strand as the gene
    @assert feature.strand == gene.gene.strand "strand error: $(string(feature)) $(string(gene.gene))"
    #check feature is within gene
    @assert feature.start >= gene.gene.start && feature.stop <= gene.gene.stop "boundary error: $(string(feature)) $(string(gene.gene))"
    if feature.type == "mRNA"
        push!(gene.mRNAs, feature)
    elseif feature.type == "tRNA"
        push!(gene.tRNAs, feature)
    elseif feature.type == "rRNA"
        push!(gene.rRNAs, feature)
    elseif feature.type == "ncRNA"
        push!(gene.ncRNAs, feature)
    elseif feature.type == "CDS"
        push!(gene.CDSs, feature)
    elseif feature.type == "exon"
        push!(gene.exons, feature)
    elseif feature.type == "intron"
        push!(gene.introns, feature)
    elseif feature.type == "misc_feature"
        push!(gene.esites, feature)
    else
        println("unknown feature type: ", feature.type)
    end
end

function add_features!(genes::Vector{Gene}, gff::DataFrame, feature_type::String)
    gff_features = filter(x->x.feature == feature_type, gff)
    for row in eachrow(gff_features)
        f = createFeature(row)
        gene = get(f.attributes, "Parent", nothing)
        if isnothing(gene)
            println("no Parent attribute for ", f.attributes["ID"])
        else
            gene = first(split(gene, "."))
            mygene = findfirst(x -> x.id == gene, genes) 
            if isnothing(mygene)
                println("no gene found matching ", f)
            else
                add_feature(genes[mygene], f)
            end
        end
    end
end

function prepact_name(position, cds_name, base, cds_position, aa, edited_aa)
    if ismissing(cds_name)
        return string(position) * "e" * string(base)
    elseif aa == AA_Gap
        return cds_name * "e" * string(base) * string(cds_position)
    else
        return cds_name * "e" * string(base) * string(cds_position) * string(aa) * string(edited_aa)
    end
end

function prepact_ref_name(species, alignments, position, cds_name, base, cds_position, aa, edited_aa)
    species_name = Dict("af" => "Af", "ar" => "Ar", "ap" => "Ap", "sm" => "Sm", "mm" => "Mm")
    if haskey(alignments,cds_name)
        ref = alignments[cds_name]["Af"]
        cds = alignments[cds_name][species_name[species]]
        aln_pos = 0
        target_pos = 0
        ref_pos = 0
        while target_pos <= cds_position
            if cds[aln_pos + 1] ≠ DNA_Gap
                target_pos += 1
            end
            aln_pos += 1
            if ref[aln_pos] ≠ DNA_Gap
                ref_pos += 1
            end
        end
        return cds_name * "e" * string(base) * string(ref_pos-1) * string(aa) * string(edited_aa) * "rAfi"
    else
        return ""
    end
end

function splice!(orfs, orf_idxs, editedseqs, gene::Gene)
    refseqs = editedseqs[gene.gene.contig]
    orfnames = unique(get.(getproperty.(gene.CDSs, :attributes), "Name", nothing))
    if gene.gene.strand == '+'
        for orfname in orfnames
            orf = LongDNA{4}()
            exons = sort(gene.CDSs[findall(x->x.attributes["Name"] == orfname, gene.CDSs)]; by = x->x.start)
            for exon in exons
                append!(orf, first(refseqs)[exon.start:exon.stop])
                idxs = haskey(orf_idxs, orfname) ? orf_idxs[orfname] : Int[]
                orf_idxs[orfname] = mod1.(append!(idxs, exon.start:exon.stop), length(first(refseqs)))
            end
            orfs[orfname] = LongRNA{4}(orf)
        end
    else
        for orfname in orfnames
            orf = LongDNA{4}()
            exons = sort(gene.CDSs[findall(x->x.attributes["Name"] == orfname, gene.CDSs)]; by = x->x.start, rev = true)
            l = length(last(refseqs))
            for exon in exons
                append!(orf, last(refseqs)[l-exon.stop+1:l-exon.start+1])
                idxs = haskey(orf_idxs, orfname) ? orf_idxs[orfname] : Int[]
                orf_idxs[orfname] = mod1.(append!(idxs, reverse(exon.start:exon.stop)), l)
            end
            orfs[orfname] = LongRNA{4}(orf)
        end
    end
end

function rc(x,l)
    l-x+1
end

function calculate_codon_position(refseqs, orf_idxs, gene::Gene, cds::Feature, position::Int)
    idxs = orf_idxs[cds.attributes["Name"]]
    cds_position = findfirst(==(position), idxs)
    if !isnothing(cds_position)
        codon_no = ceil(Int, cds_position/3.0)
        codon_position = cds_position - 3 * (codon_no - 1)
        refseq = cds.strand == '+' ? first(refseqs[cds.contig]) : last(refseqs[cds.contig])
        codon_idxs = idxs[3*codon_no-2:3*codon_no]
        if cds.strand == '-'
            codon_idxs = rc.(codon_idxs, length(refseq))
        end
        codon = LongRNA{4}(refseq[codon_idxs])
        return cds_position, codon_no, codon_position, codon
    end
    gene_position = gene.gene.strand == '+' ? position - gene.gene.start + 1 : gene.gene.stop - position + 1
    return gene_position, 0, 0, LongRNA{4}("---")
end

struct CircularSequence
    length::Int32
    sequence::LongDNA{4}
    function CircularSequence(seq::LongDNA{4})
        new(length(seq), append!(seq, LongSubSeq(seq, 1:length(seq)-1)))
    end
end

@inline Base.length(cs::CircularSequence) = cs.length
@inline Base.getindex(cs::CircularSequence, i::Integer) = @inbounds getindex(cs.sequence, mod1(i, cs.length))

function Base.getindex(cs::CircularSequence, v::Vector{<:Integer})
    bases = Vector{DNA}(undef, length(v))
    for (i,idx) in enumerate(v)
        bases[i] = cs[idx]
    end
    bases
end

function Base.getindex(cs::CircularSequence, r::UnitRange{<:Integer})
    @assert length(r) <= length(cs)
    if r.start > length(cs) || r.start < 1
        r = range(mod1(r.start, cs.length); length=length(r))
    end
    return LongSubSeq(cs.sequence, r)
end

function Base.setindex!(cs::CircularSequence, base::DNA, index::Integer)
    cs.sequence[index] = base
    cs.sequence[index + cs.length] = base
end

function BioSequences.reverse_complement(cs::CircularSequence)
    return CircularSequence(BioSequences.reverse_complement(cs.sequence[1:cs.length]))
end

function circular_in(x::Integer, interval::UnitRange{<:Integer}, circumference::Integer)
    if x ∈ interval || x + circumference ∈ interval
        return true
    else
        return false
    end
end

function main(ARGS)
    species = ARGS[1]
    organelle = ARGS[2]
    target = species * "_" * organelle
    println(target)

    refseqs = organelle == "cp" ? Dict{String, Tuple{CircularSequence, CircularSequence}}() : Dict{String, Tuple{LongDNA{4}, LongDNA{4}}}()
    FASTA.Reader(open("../data/$organelle/genomes/$target.fasta")) do infile
        refs = FASTA.Record[]
        for record in infile
            seq = FASTA.sequence(LongDNA{4}, record)
            if organelle == "cp"; seq = CircularSequence(seq); end
            refseqs[identifier(record)] = (seq, reverse_complement(seq))
        end
    end

    gff = CSV.File("../data/$organelle/genomes/$target.gff", comment = "#", header = ["sequence", "software", "feature", "start", "stop", "score", "strand", "phase", "attributes"]) |> DataFrame
    filter!(x -> x.feature ≠ "misc_feature", gff)

    gffgenes = filter(x->x.feature == "gene", gff)

    genes = Gene[]
    for row in eachrow(gffgenes)
        gene = createFeature(row)
        id = get(gene.attributes, "Name", nothing)
        if isnothing(id)
            println("no ID for ", gene.attributes["Name"])
        else
            push!(genes, Gene(gene.attributes["Name"], gene, Feature[], Feature[], Feature[], Feature[], Feature[], Feature[], Feature[], Feature[]))
        end
    end

    add_features!(genes, gff, "tRNA")
    add_features!(genes, gff, "rRNA")
    add_features!(genes, gff, "ncRNA")
    add_features!(genes, gff, "CDS")
    for gene in genes
        if gene.gene.strand == '-'
            reverse!(gene.CDSs)
        end
    end

    insites =  CSV.File("../data/$organelle/edit_sites/$target" * "_edits.tsv", comment="#") |> DataFrame

    if organelle == "cp"
        insites.refID = replace.(insites.refID, "_1kb"=>"")
        #assumes single ref seq
        refseqlength = length(first(refseqs[insites.refID[1]]))
        extension_sites = filter(x->x.pos > refseqlength, insites)
        for extensionsite in eachrow(extension_sites)
            wrapsiteview = @view(insites[insites.pos .== extensionsite.pos - refseqlength, :])
            @assert nrow(wrapsiteview) > 0 "no site equivalent to $extensionsite"
            targetsite = first(wrapsiteview)
            for col in 4:16
                targetsite[col] += extensionsite[col]
            end
            targetsite.editp = targetsite.edited/(targetsite.edited + targetsite.unedited)
            targetsite.binomial = (targetsite.binomial + extensionsite.binomial)/2.0
            targetsite.bBH = (targetsite.bBH + extensionsite.bBH)/2.0
        end
        filter!(x->x.pos <= refseqlength, insites)
    end

    editedseqs = deepcopy(refseqs)
    for site in eachrow(insites)
        seqs = editedseqs[site.refID]
        @assert DNA(site.ref[1]) == first(seqs)[site.pos] "$site\t$(first(seqs)[site.pos])"
        if site.ref == "C"
            first(seqs)[site.pos] = DNA_T
        elseif site.ref == "T"
            first(seqs)[site.pos] = DNA_C
        elseif site.ref == "A"
            last(seqs)[length(last(seqs)) - site.pos + 1] = DNA_C
        elseif site.ref == "G"
            last(seqs)[length(last(seqs)) - site.pos + 1] = DNA_T
        end
    end

    orfs = Dict{String, LongRNA{4}}()
    orf_idxs = Dict{String, Vector{Int}}()

    for g in genes
        splice!(orfs, orf_idxs, editedseqs, g)
    end

    alignments = Dict{String,Dict{String, LongDNA{4}}}()
    aln_files = filter(x -> endswith(x, ".aln"), readdir("../data/$organelle/alignments"; join = true))
    for f in aln_files
        seqs = Dict{String, LongDNA{4}}()
        FASTA.Reader(open(f)) do infile
            for r in infile
                seqs[identifier(r)] = FASTA.sequence(LongDNA{4}, r)
            end
        end
        gene = first(split(basename(f), "."))
        alignments[gene] = seqs
    end

    editing = Dict(DNA_A => RNA_G, DNA_C => RNA_U, DNA_G => RNA_A, DNA_T => RNA_C)

    esites = DataFrame(contig = String[], strand = Char[], position = Int[], id = String[], ref_id = String[], reference_base = DNA[], edited_base = RNA[], proportion_edited = Float64[], gene = String[], cds_position = Int[], codon_position = Int[],
        codon = LongRNA{4}[], edited_codon = LongRNA{4}[], aa = AminoAcid[], edited_aa = AminoAcid[], synonymous = Union{Missing, Bool}[], creates_start = Bool[], creates_stop = Bool[], removes_stop = Bool[],
        preceding_base = RNA[], subsequent_base = RNA[])
    for row in eachrow(insites)
        contig = row.refID
        position = row.pos
        refbase = DNA(row.ref[1])
        strand = '-'
        reference_base = complement(refbase)
        if refbase == DNA_C || refbase == DNA_T
            strand = '+'
            reference_base = refbase
        end
        edited_base = editing[reference_base]
        if refbase == DNA_C
            @assert edited_base == RNA_U
        elseif refbase == DNA_T
            @assert edited_base == RNA_C
        end
        proportion_edited = row.editp
        mygene = missing
        cds_name = ""
        cds_position = codon_number = codon_position = 0
        blank = LongRNA{4}("---")
        codon = blank
        for g in genes
            refseq = editedseqs[g.gene.contig]
            if row.refID == g.gene.contig && circular_in(position, g.gene.start:g.gene.stop, length(refseq)) && strand == g.gene.strand
                mygene = g
                for cds in g.CDSs
                    cds_position, codon_number, codon_position, codon = calculate_codon_position(refseqs, orf_idxs, mygene, cds, position)
                    if codon_number > 0
                        cds_name = cds.attributes["Name"]
                        break
                    end 
                end
                break
            end
        end
        ismissing(mygene) && continue
        edited_codon = LongRNA{4}([RNA_Gap, RNA_Gap, RNA_Gap])
        aa = AA_Gap
        edited_aa = AA_Gap
        synonymous = missing
        creates_start = creates_stop = removes_stop = false
        preceding_base = subsequent_base = RNA_Gap
        id = ref_id = ""
        if codon ≠ blank
            edited_codon = orfs[cds_name][3 * codon_number - 2:3 * codon_number]
            edited_codon[codon_position] = edited_base
            aa = BioSequences.translate(codon)[1]
            edited_aa = BioSequences.translate(edited_codon)[1]
            synonymous = aa == edited_aa
            creates_start = codon_number == 1 && (edited_codon ∈ [rna"AUG", rna"GUG"])
            creates_stop = cds_position == length(orfs[cds_name])-2 && (edited_codon ∈ [rna"UAA", rna"UAG", rna"UGA"])
            removes_stop = codon ∈ rna_stops
            preceding_base = orfs[cds_name][cds_position - 1]
            subsequent_base = orfs[cds_name][cds_position + 1]
            id = prepact_name(position, cds_name, edited_base, cds_position, aa, edited_aa)
            #print(id)
            ref_id = prepact_ref_name(species, alignments, position, cds_name, edited_base, cds_position, aa, edited_aa)
            #print("\t", ref_id, "\n")
        else
            preceding_base = strand == '+' ? first(refseqs[contig])[position - 1] : last(refseqs[contig])[position + 1]
            subsequent_base = strand == '+' ? first(refseqs[contig])[position + 1] : last(refseqs[contig])[position - 1]
            id = prepact_name(position, mygene.id, edited_base, cds_position, aa, edited_aa)
            #print(id)
            ref_id = ""
            #print("\t", ref_id, "\n")
        end
        push!(esites, (contig, strand, position, id, ref_id, reference_base, edited_base, proportion_edited, cds_name, cds_position, codon_position, codon, edited_codon, aa, edited_aa,
        synonymous, creates_start, creates_stop, removes_stop, preceding_base, subsequent_base))
    end
    CSV.write("../data/$organelle/edit_sites/" * target * "_" * "editing_events.tsv", esites; delim='\t')
    nothing
end
@main
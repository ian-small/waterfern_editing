using CSV, DataFrames, BioSequences, FASTX, Plots, Plots.Measures, Serialization, ProgressMeter

function map2ref(ref::Vector{DNA}, genes)
    refpos = zeros(Int, length(ref))
    count = 0
    currentgene=""
    for (n,(g,base)) in enumerate(zip(genes,ref))
        if g ≠ currentgene
            count = 0
            currentgene = g
        end
        if base ≠ DNA_Gap
            count += 1
            refpos[n] = count
        end
    end
    refpos
end

function onlyY(bases)
    if all(x->x ∈ [DNA_C, DNA_T], bases)
        return true
    else
        return false
    end
end

function Yissynonymous(genes::Vector{AbstractString}, bases::Vector{DNA})
    syn = falses(length(bases))
    codonbases = DNA[]
    gene = ""
    for (i,(g, b)) in enumerate(zip(genes, bases))
        if g ≠ gene
            empty!(codonbases)
            gene = g
        end
        if b ≠ DNA_Gap
            push!(codonbases, b)
        end
        if length(codonbases) == 3
            codon = LongDNA{4}(codonbases)
            aa = BioSequences.translate(codon)[1]
            for codonposition in 1:3
                codon = LongDNA{4}(codonbases)
                issyn = true
                codon[codonposition] = DNA_C
                if BioSequences.translate(codon)[1] ≠ aa
                    issyn = false
                end
                codon[codonposition] = DNA_T
                if BioSequences.translate(codon)[1] ≠ aa
                    issyn = false
                end
                syn[i-3+codonposition] = issyn
            end
            empty!(codonbases)
        end
    end
    syn
end

function parse_uid(uid::AbstractString)
    m = match(r"([A-Za-z0-9]+)e([C|U])([0-9]+)", uid)
    gene = m.captures[1]
    type = first(m.captures[2])
    refpos = parse(Int, m.captures[3])
    gene, type, refpos
end

function differs_by_mutation(edittype::Char, dna::Vector{DNA}, editing::BitVector)
    for (d, e) in zip(dna, editing)
        if edittype == 'U'
            if (d == DNA_C && ~e)
                return false
            end
        elseif (d == DNA_T && ~e)
            return false
        end
    end
    true
end

# Trajectories are based on this tree topology, branches indexed 1-8 from left to right:
#   Af Ar Ap Sm Mm
#   |__|  |  |  |
#     |___|  |  |
#       |____|  |
#          |____|
#             |

const branch_lengths = Dict{Int, Float64}(7=>64.0,8=>153.0,5=>38.3,6=>89.0,1=>3.1,2=>3.1,3=>47.6,4=>50.7, 9=>0.0)

function trajectory2pattern(trajectory::BitVector)
    pattern = fill(trajectory[9], 5)
    pattern[5] = pattern[5] ⊻ trajectory[8]
    pattern[4] = pattern[4] ⊻ trajectory[7] ⊻ trajectory[6]
    pattern[3] = pattern[3] ⊻ trajectory[7] ⊻ trajectory[5] ⊻ trajectory[4]
    pattern[2] = pattern[2] ⊻ trajectory[7] ⊻ trajectory[5] ⊻ trajectory[3] ⊻ trajectory[2]
    pattern[1] = pattern[1] ⊻ trajectory[7] ⊻ trajectory[5] ⊻ trajectory[3] ⊻ trajectory[1]
    pattern
end

function extend(perms::Vector{BitVector})
    newperms = Vector{BitVector}()
    for perm in perms
        push!(newperms, push!(copy(perm), 0))
        push!(newperms, push!(perm, 1))
    end
    newperms
end

function bitperm(length::Int)
    perms = Vector{BitVector}()
    push!(perms, BitVector())
    for n in 1:length
        perms = extend(perms)
    end
    perms
end

#ML approach
#likelihood of transition event occurring on a particular branch is the rate * branch length
#likelihood of a specific trajectory is the product of the likelihoods of the sequential events
#likelihood of a specific pattern is the sum of the likelihoods of the trajectories leading to that pattern
#likelihood of overall data (set of observed patterns) is the product of the likelihoods of each trajectory

function trajectory_likelihood(ancestral::Float64, trajectory::BitVector, TtoC::Float64, CtoT::Float64)

    function branch_likelihood(lca::Bool, branch::Int)
        likelihood = 1.0
        time = branch_lengths[branch]
        if lca
            if trajectory[branch] 
                likelihood = 1.0 - (ℯ^-(CtoT * time))
            else
                likelihood = ℯ^-(CtoT * time)
            end
        else
            if trajectory[branch]
                likelihood = 1.0 - (ℯ^-(TtoC * time))
            else                            
                likelihood = ℯ^-(TtoC * time)
            end
        end
        likelihood
    end

    likelihood = trajectory[9] ? ancestral : 1 - ancestral
    lca = trajectory[9]
    likelihood *= branch_likelihood(lca, 8)
    likelihood *= branch_likelihood(lca, 7)
    lca = trajectory[9] ⊻ trajectory[7]
    likelihood *= branch_likelihood(lca, 6)
    likelihood *= branch_likelihood(lca, 5)
    lca = trajectory[9] ⊻ trajectory[7] ⊻ trajectory[5]
    likelihood *= branch_likelihood(lca, 4)
    likelihood *= branch_likelihood(lca, 3)
    lca = trajectory[9] ⊻ trajectory[7] ⊻ trajectory[5] ⊻ trajectory[3]
    likelihood *= branch_likelihood(lca, 2)
    likelihood *= branch_likelihood(lca, 1)
    likelihood
end

# first partition is default to which default TtoC, CtoT rates apply; all subsequent partitions have specific rate
struct PartitionData
    partitions::Vector{Vector{Int}}     # each partition is a vector of indices into the editing patterns
    partition_directions::BitVector     # 1 indicates partition-specific rate is TtoC, 0 indicates it is CtoT
end

# rates: generic TtoC, generic CtoT, specific rate for each partition
function loglikelihood_f!(rates::Vector{Float64}, data::PartitionData, ancestral::Float64,
            trajectories::Vector{BitVector}, trajectory_likelihoods::Matrix{Float64},
            trajectoriesbypattern::Vector{Vector{Int}}, pattern_likelihoods::Matrix{Float64}
            )::Float64
    
    for p in 1:length(data.partitions)
        TtoC = rates[1]
        CtoT = rates[2]
        if p > 1
            if data.partition_directions[p-1]
                TtoC = rates[p+1]
            else
                CtoT = rates[p+1]
            end
        end
        for (i,t) in enumerate(trajectories)
            trajectory_likelihoods[p,i] = max(0.0, trajectory_likelihood(ancestral, t, TtoC, CtoT))
        end
        for (i, t_idxs) in enumerate(trajectoriesbypattern)
            pl = 0.0
            for idx in t_idxs
                pl += trajectory_likelihoods[p,idx]
            end
            pattern_likelihoods[p,i] = pl
        end
    end
    loglikelihood = 0.0
    for p in 1:length(data.partitions)
        partition = data.partitions[p]
        likelihoods = pattern_likelihoods[p,:]
        loglikelihood += sum(log.(getindex(likelihoods, partition)))
    end
    loglikelihood
end

function profile_likelihoods(pdata::PartitionData, searchspaces, resolution::Int, ancestralrange,
    trajectories::Vector{BitVector}, trajectoriesbypattern::Vector{Vector{Int}})
    num_partitions = length(pdata.partitions)
    num_rates = num_partitions + 1
    profiles = Array{Float64}(undef, fill(resolution, num_rates + 1)...)
    @showprogress Threads.@threads for i in eachindex(profiles)
        ancestral = 0.5 #default, isn't used
        rates = Vector{Float64}(undef, num_rates)
        trajectory_likelihoods = Matrix{Float64}(undef, num_partitions, length(trajectories))
        pattern_likelihoods = Matrix{Float64}(undef, num_partitions, length(trajectoriesbypattern))
        for (s, j) in enumerate(Tuple(CartesianIndices(profiles)[i]))
            if s <= length(searchspaces)
                rates[s] = searchspaces[s][j]
            else
                ancestral = ancestralrange[j]
            end
        end
        profiles[i] = loglikelihood_f!(rates, pdata, ancestral,
            trajectories, trajectory_likelihoods,
            trajectoriesbypattern, pattern_likelihoods)
    end
    profiles
end

function main(ARGS)
    genes = AbstractString[]
    species = ["Af", "Ap", "Ar", "Sm", "Mm"]
    organelle = "cp" # change to "mt" for mitochondrial profiles
    sequences = Dict{String, LongDNA{4}}()
    for s in species
        sequences[s] = LongDNA{4}()
    end
    geneoffsets = [0]

    for filename in sort(filter(x->endswith(x, ".aln"), readdir("data/$organelle/alignments"; join = true)))
        gene = first(split(basename(filename), "."))
        FASTA.Reader(open(filename)) do alnfile
            for (n, record) in enumerate(alnfile)
                seq = sequence(LongDNA{4}, record)
                push!(geneoffsets, last(geneoffsets)+length(seq))
                if n == 1
                    append!(genes, fill(gene, length(seq)))
                end
                sequences[identifier(record)] = append!(sequences[identifier(record)], seq)
            end
        end
    end

    df = DataFrame(pos=collect(1:length(genes)),gene=genes,Af=collect(sequences["Af"]),Ap=collect(sequences["Ap"]),Ar=collect(sequences["Ar"]),Sm=collect(sequences["Sm"]),Mm=collect(sequences["Mm"]))
    df.refpos = map2ref(df.Af, genes)
    df.onlyY = onlyY.(eachrow(df[:,3:7]))
    filter!(x->x.onlyY, df)
    println("$(nrow(df)) DNA sites considered")

    synAf = Yissynonymous(df.gene, df.Af)
    synAp = Yissynonymous(df.gene, df.Ap)
    synAr = Yissynonymous(df.gene, df.Ar)
    synSm = Yissynonymous(df.gene, df.Sm)
    synMm = Yissynonymous(df.gene, df.Mm)

    df.Yissynonymous = synAf .& synAp .& synAr .& synSm .& synMm

    #add editing sites
    df.editsite = falses(nrow(df))
    df.edittype = fill('.', nrow(df))
    df.creates_start = falses(nrow(df))
    df.creates_stop = falses(nrow(df))
    df.removes_stop = falses(nrow(df))
    df.editing_differs_by_mutation = trues(nrow(df))
    editing_sites = CSV.File("data/$organelle/edit_sites/all_$organelle" * "_sites.tsv") |> DataFrame
    for site in eachrow(editing_sites)
        gene, type, refpos = parse_uid(site.uid)
        if nrow(df[(df.gene .== gene) .& (df.refpos .== refpos), :]) == 0
            continue
        end
        df[(df.gene .== gene) .& (df.refpos .== refpos), :editsite] .= true
        df[(df.gene .== gene) .& (df.refpos .== refpos), :edittype] .= type
        dna = first(df[(df.gene .== gene) .& (df.refpos .== refpos), :])
        df[(df.gene .== gene) .& (df.refpos .== refpos), :editing_differs_by_mutation] .= differs_by_mutation(type, [dna.Af,dna.Ar,dna.Ap,dna.Sm,dna.Mm], .~ismissing.([site.af,site.ar,site.ap,site.sm,site.mm]))
        df[(df.gene .== gene) .& (df.refpos .== refpos), :creates_start] .= site.creates_start
        df[(df.gene .== gene) .& (df.refpos .== refpos), :creates_stop] .= site.creates_stop
        df[(df.gene .== gene) .& (df.refpos .== refpos), :removes_stop] .= site.removes_stop
    end
    filter!(x->x.editsite || x.Yissynonymous, df)
    filter!(x->x.editing_differs_by_mutation, df)

    println("$(nrow(df)) sites retained, of which $(sum(df.editsite)) sites are edited")

    println("constructing evolutionary trajectories...")
    trajectories = bitperm(9)
    patterns = bitperm(5)
    trajectoriesbypattern = [Int[] for i in 1:length(patterns)]
    for (i,t) in enumerate(trajectories)
        p_idx = findfirst(==(trajectory2pattern(t)), patterns)
        push!(trajectoriesbypattern[p_idx], i)
    end

    # generate editing patterns and store tham as the index of the matching pattern in patterns
    println("collating editing patterns across species...")
    epatterns = Int[]
    for site in eachrow(df)
        p = ==(DNA_C).([site.Af,site.Ar,site.Ap,site.Sm,site.Mm])
        push!(epatterns, findfirst(==(p), patterns))
    end
    df.epattern = epatterns

    #unpartitioned
#=  pdata=PartitionData([df.epattern],BitVector[])
    resolution = 100
    searchspace = logrange(1e-4, 1e-2, resolution)
    searchspaces = fill(searchspace, length(pdata.partitions) + 1)
    profiles = profile_likelihoods(pdata, searchspaces, resolution,
        ancestral, trajectories, trajectoriesbypattern)
    serialize("unpartitioned.bin", profiles)

    #partitioned
    pdata = PartitionData([df[df.edittype .== '.', :epattern], df[df.edittype .== 'U', :epattern], df[df.edittype .== 'C', :epattern]],BitVector([0,1]))
    @assert sum(length.(pdata.partitions)) == nrow(df)
    resolution = 10
    searchspace = logrange(1e-4, 1e-2, resolution)
    searchspaces = fill(searchspace, length(pdata.partitions) + 1)
    profiles = profile_likelihoods(pdata, searchspaces, resolution,
        ancestral, trajectories, trajectoriesbypattern)
    serialize("partitioned.bin", profiles) =#

    #fully partitioned
    pdata = PartitionData(
    [df[df.edittype .== '.', :epattern],
    df[(df.edittype .== 'U') .& .~df.Yissynonymous .& .~df.creates_start .& .~df.creates_stop, :epattern],
    df[df.creates_start, :epattern],
    df[df.creates_stop, :epattern],
    df[(df.edittype .== 'C') .& .~df.Yissynonymous .& .~df.removes_stop, :epattern],
    df[df.removes_stop, :epattern]],
    BitVector([0,0,0,0,1,1,1]))

    if organelle == "cp"
        resolution = 11
        ancestral = range(0.39,0.44; length = resolution)
        searchspaces = Base.LogRange[]
        push!(searchspaces, logrange(10^-3.4, 10^-3.2, resolution))
        push!(searchspaces, logrange(10^-3.3, 10^-3, resolution))
        push!(searchspaces, logrange(10^-2.4, 10^-2.2, resolution))
        push!(searchspaces, logrange(10^-2.8, 10^-2.3, resolution))
        push!(searchspaces, logrange(10^-2.9, 10^-1.8, resolution))
        push!(searchspaces, logrange(10^-3.2, 10^-1.8, resolution))
        push!(searchspaces, logrange(10^-3.6, 10^-2.66, resolution))
        profiles = profile_likelihoods(pdata, searchspaces, resolution,
            ancestral, trajectories, trajectoriesbypattern)
        serialize("fully_partitioned_cp_zoom.bin", profiles)
    elseif organelle == "mt"
        resolution = 11
        ancestral = range(0.40,0.44; length = resolution)
        searchspaces = Base.LogRange[]
        push!(searchspaces, logrange(10^-3.7, 10^-3.45, resolution))
        push!(searchspaces, logrange(10^-5.0, 10^-4.2, resolution))
        push!(searchspaces, logrange(10^-3.2, 10^-2.9, resolution))
        push!(searchspaces, logrange(10^-6.0, 10^-2.95, resolution))
        push!(searchspaces, logrange(10^-6.0, 10^-3.0, resolution))
        push!(searchspaces, logrange(10^-2.6, 10^-2.0, resolution))
        push!(searchspaces, logrange(10^-3.3, 10^-2.9, resolution))
        profiles = profile_likelihoods(pdata, searchspaces, resolution,
            ancestral, trajectories, trajectoriesbypattern)
        serialize("fully_partitioned_mt_zoom.bin", profiles)
    end
end
@main



    
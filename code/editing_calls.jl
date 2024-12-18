#read in packages to be used
using CSV
using DataFrames
using FASTX
using HypothesisTests
using MultipleTesting

function main(counts_file::String, output_file::String)

    #the CSV package is creating a Table from the .count file and the DataFrame package is creating a DataFrame from the Table. 
    rcounts = DataFrame(CSV.File(counts_file, delim="\t", comment ="#"))

    #drop the :total column from the DataFrame counts
    rcounts = select!(rcounts, Not(:total))

    genomelength = nrow(rcounts)

    #forward count total
    rcounts.f_total = rcounts.fA .+ rcounts.fC .+ rcounts.fG .+ rcounts.fU

    #reverse count total
    rcounts.r_total = rcounts.rA .+ rcounts.rC .+ rcounts.rG .+ rcounts.rU

    #create extra columns with zeros for edited and unedited counts shown on both the forward and reverse strands.
    rcounts.edited = zeros(Int, nrow(rcounts))
    rcounts.unedited = zeros(Int, nrow(rcounts))

    #forward strand C-U editing
    #total C edited on forward strand (C on reference) is equal to the total U on the forward strand (when ref = C).
    rcounts[rcounts.ref .== "C", :edited] = rcounts[rcounts.ref .== "C", :fU]
    #total C unedited on the forward strand (C on reference) is equal to the total on the forward strand, minus the number of edited shown on the forward strand (when ref = C).
    rcounts[rcounts.ref .== "C", :unedited] = rcounts[rcounts.ref .== "C", :f_total] .- rcounts[rcounts.ref .== "C", :edited]

    #forward strand U-C editing
    #total U edited on forward strand (T on reference) is equal to the total C on the forward strand (when ref = T).
    rcounts[rcounts.ref .== "T", :edited] = rcounts[rcounts.ref .== "T", :fC]
    #total U unedited on the forward strand (T on reference) is equal to the total on the forward strand, minus the number of edited shown on the forward strand (when ref = T).
    rcounts[rcounts.ref .== "T", :unedited] = rcounts[rcounts.ref .== "T", :f_total] .- rcounts[rcounts.ref .== "T", :edited]
 
    #reverse strand C-U editing
    #total C edited on the reverse strand (G on reference) is equal to the total U on the reverse strand (when ref = G). 
    rcounts[rcounts.ref .== "G", :edited] = rcounts[rcounts.ref .== "G", :rU]
    #total C unedited on the reverse strand (G on reference) is equal to the total on the reverse strand, minus the number of edited shown on the reverse strand (when ref = G).
    rcounts[rcounts.ref .== "G", :unedited] = rcounts[rcounts.ref .== "G", :r_total] .- rcounts[rcounts.ref .== "G", :edited]
   
    #reverse strand U-C editing
    #total U edited on the reverse strand (A on reference) is equal to the total C on the reverse strand (when ref = A). 
    rcounts[rcounts.ref .== "A", :edited] = rcounts[rcounts.ref .== "A", :rC]
    #total U unedited on the reverse strand (A on reference) is equal to the total on the reverse strand, minus the number of edited shown on the reverse strand (when ref = A).
    rcounts[rcounts.ref .== "A", :unedited] = rcounts[rcounts.ref .== "A", :r_total] .- rcounts[rcounts.ref .== "A", :edited]


    #total counts = the sum of edited and unedited shown on the forward strand.
    rcounts.total = rcounts.edited .+ rcounts.unedited
    #what proportion of the total is edited (shown on the forward strand).
    rcounts.editp = rcounts.edited ./ rcounts.total

    #run a binomial test with a pvalue of 0.05 for Ph1 counts, single-tailed to the right.
    rcounts.binomial = pvalue.(BinomialTest.(rcounts.edited, rcounts.total, 0.05); tail = :right)
    #use BenjaminiHochberg to adjust the individual p-value for each gene to keep the overall error rate (or false positive rate) to less than or equal to the user-specified p-value cutoff or error rate.
    rcounts.bBH = adjust(rcounts.binomial, BenjaminiHochberg())

    #All editing sites after adjusted (Benjamini-Hochberg), p-value less than 0.05.
    editing_sites_binomial = rcounts[(rcounts.bBH .< 0.05),:]

    CSV.write(output_file, editing_sites_binomial, delim="\t")

end

main(ARGS...)

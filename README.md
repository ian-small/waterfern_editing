## Repository of data and code for the manuscript 'High conservation of translationally significant RNA editing sites in hyper-editing ferns implies they are not selectively neutral' by van der Guizen et al.

The data in this repository is from five water ferns from the *Salviniales* family: *Azolla filiculoides*, *Azolla rubra*, *Azolla pinnata*, *Salvinia molesta* and *Marsilea mutica*. The raw sequence reads, genome assemblies and transcriptome assemblies are not included in this repository but are available from the [SRA](https://www.ncbi.nlm.nih.gov/sra) or [GenBank](https://www.ncbi.nlm.nih.gov/genbank/). The relevant accessions are listed in `data/accessions.tsv`.

### Read mapping

RNA reads were trimmed using BBDuk from the [BBtools](https://sourceforge.net/projects/bbmap/) suite:  
`bbduk.sh in1=read1.fq.gz in2=read2.fq.gz out1=read1.trimmed.fq.gz out2=read2.trimmed.fq.gz ktrim=r k=23 mink=11 hdist=1 ftm=5 tpe tbo`  
then merged:  
`bbmerge.sh in1=read1.trimmed.fq.gz in2=read2.trimmed.fq.gz out=merged.fq.gz outu=unmerged.fq.gz qtrim=2 trimq=10,15,20, minq=12`

Trimmed reads were mapped to the target genome with bbwrap/bbmap:  
`bbwrap.sh in=merged.fq.gz, unmerged.fq.gz ref=target.fa out=target.bam mappedonly ambiguous=random`

### Calling editing sites

Counts of the RNA nucleotides present at every position in the `.bam` files were generated for each species using Pyrimid (source code in `code/pyrimid`):  
`julia code/pyrimid/src/pyrimid.jl -m 0 -u -o target.counts target.fa target.bam`  
`-m 0` allows zero non-editing mismatches (any read with mismatches to the reference that cannot be explained by editing will be discarded)  
`-u` allows for U-to-C edits to be considered as well as C-to-U edits

The resulting counts files are included in this repository in `data/cp/counts` and `data/mt/counts`.

Editing sites were identified as sites with a statistically significant probability of > 5% editing using a binomial test with Benjamini-Hochberg multiple testing correction:  
`julia code/pyrimid/editing_calls.jl target.counts target_edits.tsv`

The resulting tables are included in this repository in `data/cp/edit_sites` and `data/mt/edit_sites`.

### Analysing editing sites across species

Editing sites were assigned IDs of the form `[gene][e][edited base][position][translation from unedited codon][translation from edited codon]`, e.g. `matKeU545PL` for the C-to-U event at the 545th nucleotide of the *matK* CDS in *Azolla pinnata* that results in a proline -> leucine difference in the translated protein. They were then assigned universal IDs of the same form using the *Azolla filiculoides* sequences as references, e.g. `matKeU545PL` becomes `matKeU536PLrAfi` as nucleotide 545 of the *Azolla pinnata* *matK* CDS aligns with the nucleotide 536 of the *Azolla filiculoides* *matK* CDS. In the process, translationally significant editing sites are idebtified and flagged. The code that does these mappings is `code/gff2esites.jl` whch can be run from within Julia as follows:

`julia> include("code/gff2esites.jl")`  
`julia> for species in ["af", "ar", "ap", "sm", "mm"], organelle in ["cp", "mt"]; main([species, organelle]); end`

The tables the code generates are included in this repository as `data/cp/edit_sites/*_cp_editing_events.tsv` and `data/mt/edit_sites/*_mt_editing_events.tsv`. These tables were merged in to a single summary table for each organelle (`data/cp/edit_sites/all_cp_sites.tsv`, `data/mt/edit_sites/all_mt_sites.tsv`) using the notebook `code/merge_species.ipynb`. These merged tables were used as the raw data for Fig. 3, Fig. 4, Fig 5C, Fig. 6, Table S2, Table S3 and Table S4.

`code/site_distribution.ipynb` plots the distribution of editing sites of different classes within the first 1000 nt of organellar coding seqences. This plot is shown in Fig. S7.

### Maximum parsimony

`code/parsimony.ipynb` estimates editing site gains and losses for each organelle by the method of maximum  parsimony. Where different evolutionary scenarios are equally parsimonious, the events are distributed randomly amongst the branches in proportion to the branch lengths (i.e. longer branches are attributed more events). One such estimate is shown in Fig 3C.

### Estimation of C→T and T→C transition rates

C→T and T→C transition rates were estimated using a non-time-reversible model by a maximum likelihood approach. The parameters of the model comprised C→T and T→C transition rates for each site class (unedited synonymous T→C, unedited synonymous C→T, non-synonymous C→T, start codon creation C→T, stop codon creation C→T, non-synonymous T→C, stop codon removal T→C) and the probability that the ancestral state of the site was C. To gain an idea of the confidence we can have in the parameter estimates, we carried out a profile likelihood analysis, by varying each parameter across a range either side of the optimal value. The resulting profile likelihoods are shown in Fig. S6, and the estimates of the C→T and T→C transition rates and their confidence intervals are shown in Fig. 5A and 5B. The code for calculating the profile likelihoods is included as `code/profile_likelihoods.jl` and can be run as follows:  
`julia -t 16 code/profile_likelihoods.jl cp`  or `julia -t 16 code/profile_likelihoods.jl mt`
It is recommended to run the code using multiple threads as it requires a very large number of calculations (for 7 rate parameters plus the ancestral state parameter and a grid resolution of 11, the likelihood of the observed editing patterns must be calculated 8^11 times, i.e. ~8.6 billion times). This will take many hours, even with multi-threading.

If you are just interested in the less complex models (i.e bipartitioned instead of fully partitioned), you can comment out `partioning_scheme = "fully_partitioned"` (line 476 in `profile_likelihoods.jl`) and uncomment the lines that activate the partioning schemes you are interested in. Bipartitions will run much more quickly (a few minutes).

The data output by `profile_likelihoods.jl` is not included in the repository because it is rather voluminous (~1.7 GB), but can be generated by running the code as shown above. A [Jupyter](https://jupyter.org) notebook (using Julia code) is provided (`profile_likelihoods.ipynb`) to visualise the data. To use Julia in a Jupyter notebook you must have [IJulia](https://github.com/JuliaLang/IJulia.jl) installed. The code in the notebook was used to generate Fig. S6 and Fig. 5A and 5B.

Data output from the bipartitions can be visualised with `code/profile_likelihoods_site_positions.ipynb` and `code/profile_likelihoods_triplets.ipynb`. The calculations from this code are shown in Table S5.



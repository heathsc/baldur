# baldur
Variant (SNV and small indel) calling from ONT sequencing data on the mitochondria 

To compile you will need an up-to-date copy of rust.  This can be
installed locally following the instructions [here](https://www.rust-lang.org/learn/get-started).  
Note that if you have rust already installed you should update it
using ``rustup update`` before trying to compile ont_demult.

Clone the baldur repository and then from the ont_demult directory
use cargo to compile the application:

    cargo build --release
	 
After successful the executable will be found in target/release/.  It
should be copied somewhere where it can be found by the shell.

Once installed, basic help can be found by invoking baldur with
the -h flag.

---------------
Basic operation
---------------

baldur in invoked with a reference fasta file for the mitochondria,
and a BAM alignment file sorted in read order
(i.e., as it comes from the mapper), but there are many options to
allow specification of thresholds for base and mapping qualities,
lists of rs ids, blacklists etc.  A typical command line would be:

``baldur --rs-list bed_chr_MT.bed.gz -q20 -Q10 -M30 -P3 -T chrM.fa -o output input.bam``

The above line would run baldur using `chrM.fa` as reference, readinf from `input.bam` and 
generating an output file `output.vcf.gz`.  The other options set the MAPQ theshold to 20, the base quality threshold to 10, the maximum base quality to 30 (any quality above this will be set to 30), and the
minimum size of a homopolymer to generate a warning filter to 3

-------
Changes
-------

1.1.0 Generate large deletions calls separately from smaller deletions.  For large deletions we don't explicitly consider the non-deleted allele
and just count the number of deletions against the other alleles.  Large deletions with similar endpoints are merged together and the collection of observed
deletion alleles at each location are used to estimate empirical confidence intervals for the start and length of the allele. After calling the deletion the observations
contributing to the deletion are removed and variant calling is performed using the remaining observations.  This allows calling of variants that fall 
in the non-deleted allele.


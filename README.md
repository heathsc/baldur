# baldur
Variant (SNV and small indel) calling from ONT sequencing data on the mitochondria 

To compile you will need an up to date copy of rust.  This can be
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
lists of rs ids, blacklists etc.

``baldur --blacklist blacklist.txt --rs-list bed_chr_MT.bed.gz -M 30 -P 3 -T chrM.fa -o output input.bam``





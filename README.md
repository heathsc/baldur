## baldur
Variant (SNV, small indel and large deletion) calling from ONT sequencing data on the mitochondria 

 - [Introduction](#intro)
 - [Installation](#install)
 - [General usage](#usage)
   - [Command line options](#cli)
 - [Overview of operation](#overview)
 - [Changes](#changes)

## <a name="intro"></a>Introduction

Baldur is a utility for variant calling of single nucleotide, short indel and large deletions from ONT sequencing data on the mitochondria. 
While multiple utilities are available for variant calling from ONT data, and there is at least one pipeline specialized for calling sequence
variants from mitochondrial data, to our knowledge baldur is the first utility specialized for calling variants from ONT data on mitochondria. 

The focus of the analysis is not on calling genotypes since the mitochondria can not be considered haploid or diploid but instead *n*-ploid where *n* can be large (roughly corresponding to the 
number of mitochondria in a cell).  Instead the focus is on finding loci that are
*heteroplasmic*, and providing estimates of the heteroplasmy frequencies, along with 
confidence bounds. At any position multiple variants may be present, and there is no *a priori* assumption 
that variants at different loci form haplotypes as different variants can arise independently in different mitochondria.

The relatively high error rate of ONT sequence data (compared to Illumina sequencing) places limits on the 
detection of low frequency variants; this is particularly true for small (1-3bp) indels where the error rate can 
be from 5-10%, however with enough coverage and more recent ONT data it is possible to call SNVs in the 
range 0.5% - 1% with confidence. 

The main output of baldur is a VCF file with the called variants.  In addition, a file showing the sequencing depth and coverage across the mitochondria is produced, which is useful for
checking the sequencing output and looking for large (multi-kb) deletions.  An optional extra output file can
also be produced with the full pileup; this is a text file with 1 row per read and 1 column per base position, making
it simple to explore combinations of multiple variant loci.

## <a name="install"></a>Installation

To compile you will need an up-to-date copy of rust.  This can be
installed locally following the instructions [here](https://www.rust-lang.org/learn/get-started).  
Note that if you have rust already installed you should update it
using ``rustup update`` before trying to compile baldur.

You will also need to have [htslib](https://github.com/samtools/htslib) installed in a place 
where the rust compiler can find it.  Note that *baldur* has been tested with versions of htslib from 1.12 to 1.15, but other recent versions should also work.

Clone the repository and then from the baldur directory
use cargo to compile the application:
```
    git clone https://github.com/heathsc/baldur.git
    cd baldur
    cargo build --release
```
If you have an error during linkage about not being able to find libhts then you will need to specify the installation location of libhts
during the build process:

    RUSTFLAGS="-L /opt/share/htslib/1.14/lib/"  cargo build --release

After successful the executable will be found in target/release/.  It
should be copied somewhere where it can be found by the shell.  The file `share/blacklist.txt` should also be copied
to a generally accessible location if you intend to use the [blacklist](#blacklist) option

Once installed, basic help can be found by invoking baldur with
the -h flag.

## <a name="usage"></a>General usage

Baldur is invoked with a reference fasta file for the mitochondria,
and a BAM alignment file sorted in read order
(i.e., as it comes from the mapper), but there are many options to
allow specification of thresholds for base and mapping qualities,
lists of rs ids, blacklists etc. A minimal command line invocation would therefore be:
```
baldur -T chrM.fa.gz in.bam
```
And a more typical command line showing multiple options would be:
```
baldur --rs-list bed_chr_MT.bed.gz --blacklist share/backlist.txt -q20 -Q10 -M30 -T chrM.fa -o output input.bam
```

The above line would run baldur using `chrM.fa` as reference, reading from `input.bam` and 
generating an output file `output.vcf.gz`.  The other options set the MAPQ theshold to 20, 
the base quality threshold to 10, the maximum base quality to 30 (any quality above this will
be set to 30), and specify that SNP rs ids can be found from the file `bed_chr_MT.bed.gz` and that a blacklist
of positions to be filtered should be read from the file `share/blacklist.txt`.

The circular nature of the mitochondria can make read mapping challenging, and several strategies can be used
to get around this.  One of the more common is to use a *rotated reference* where the origin of the reference
has been shifted by a certain number of bases (which can help map reads that span the original origin).  

Baldur has two behaviours that help the handling of alignments from such shifted references.
1. There is an option (`--adjust`) to internally adjust the map positions of the alignments in the input BAM file (to compensate for the rotation).
2. Internally baldur treats each position as modulo the genome length so a position beyond the end of the genome wraps around to the beginning.

### <a name="cli"></a>Command line options

Baldur has many command line options for controlling the operation of the calling process.

| Short                    | Long                 | Description                                                 | Default        |
|--------------------------|----------------------|-------------------------------------------------------------|----------------|
| q                        | mapq-threshold       | MAPQ threshold                                              | 0              |
| Q                        | qual-threshold       | Minimum base quality                                        | 0              |
| M                        | max-qual             | Quality values for bases are capped at this value           | 30             |
| I                        | max-indel-qual       | Quality values for indels are capped at this value          | 20             |
| P                        | homopolymer-limit    | Minimum size of homopolymer runs flagged as problematic     | 4              |
|                          | snv-thresholds       | Hard and soft limits for SNV reporting                      | 0.0005, 0.0025 |
|                          | indel-thresholds     | Hard and soft limits for indel reporting                    | 0.025, 0.1     |
|                          | small-deletion-limit | Maximum size for a small (explicit) deletion                | 64             |
|||||
| a                        | adjust               | Adjustment to genomic position of alignments                | 0              |
| r                        | region               | Genomic region to consider                                  |                |
|||||
| o                        | output-prefix        | Prefix for output files                                     | baldur         |
| n                        | sample               | Sample name (in VCF file)                                   | SAMPLE         |
| T                        | reference            | Reference FASTA file (**REQUIRED**)                         |                |
| l                        | loglevel             | Set log level (none, error, warn, info, debug, trace)       | info           |
| <a name="blacklist"></a> | blacklist            | BED file with list of blacklisted sites                     |                |
|                          | rs-list              | BED file with dbSNP rs identifiers                          |                |
|                          | view          | Generate pileup view file|                |
|                          | no-call              | Do not perform variant calling                              |                |

## <a name="overview"></a>Overview of workflow

- Read in alignments filtering on region and MAPQ
- Generate internal pileup representation of all reads
- Output pileup view if required
- Output sequencing depth report
- Estimate single site allele frequencies (SNVs and indels)
- Identify long deletions
- Merge consecutive short indels and estimate multisite allele frequencies
- Generate VCF

## <a name="changes"></a>Changes

- 1.1.7 Switch to r_htslib 0.4
- 1.1.6 Switch to r_htslib 0.3
- 1.1.5 Switch compress_io version to 0.5.
- 1.1.5 Change documentation to suggest using https rather than ssh
- 1.1.4 Switch to using compress_io from crates.io
- 1.1.3 Add allele count output for long deletion alleles.  Correct bug where low quality deletions were being counted as non-deleted alleles.
- 1.1.2 Fix bug in merging of deletion alleles (introduced in 1.1.0)
- 1.1.1 Change default quality thresholds for bases and indels to 30 and 20 respectively, and change default indel threshold to 0.025, 0.1 (hard, soft).
- 1.1.0 Generate large deletions calls separately from smaller deletions.  For large deletions we don't explicitly consider the non-deleted allele
and just count the number of deletions against the other alleles.  Large deletions with similar endpoints are merged together and the collection of observed
deletion alleles at each location are used to estimate empirical confidence intervals for the start and length of the allele. After calling the deletion the observations
contributing to the deletion are removed and variant calling is performed using the remaining observations.  This allows calling of variants that fall 
in the non-deleted allele.


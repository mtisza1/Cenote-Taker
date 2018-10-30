### _Cenote-Taker Overview_

***** Before attempting to use Cenote-Taker, know that this is a work in progress. While it works great for me, it will likely be very difficult to port it to your server based on the code as is. If you are in dire need of Cenote-Taker, please contact Mike. *****

Few of the complete viral genomes assembled from metagenomic deep sequencing reads have been deposited on NCBI's GenBank database, contributing to a lack of viral reference genomes available to the scientific community for analysis and comparison.  The primary reasons for not depositing genomes include challenges in gene annotation, dealing with file format conversions, and determining which contigs represent complete genomes. Relatedly, the arduousness in manually curating maps of Open Reading Frames (ORFs) in novel viral genomes and contigs hinders analysis.

Cenote-Taker (a portmanteau of _Cenote_, a naturally occuring circular water pool, and _note-taker_) was designed to overcome the most challenging obstacles in high throughput viral (and plasmid) genome/sequence annotation and facilitation of GenBank deposition of annotated genomes. Cenote-Taker should work for annotation of genomes closely related to extant reference genomes, but is specifically designed for annotation of sequences that are highly divergent from known sequences.

For users interested in discovering which sequences amongst their metagenomic contigs are circular, the '-default' option can take contigs directly from de novo assembly (ideally assembled with SPADes with the '--plasmid' option), scan each sequence for terminally redundant ends (indicating a complete circular DNA molecule), and annotate all ORFs on circular sequences.

For users who have a set of contigs that they have already determined circularity for and 'wrapped', the options '-given_circular' or '-needs_rotation' can be submitted, depending on whether the user has rotated the contig to exclude ORFs from crossing the 'wrap point' of the genomes.

For users who have a set of linear genomes or contigs, the '-given_linear' option can be used to annotate their sequences.

GenBank requires that genome submissions include metadata about each genome. Therefore, Cenote-taker gives users the opportunity to include that information to be incorporated into the output.
 

### _Cenote-Taker's processes_

### **Part 1: Contig processing**
_-default option_
(1) The specified .fasta file will be parsed to remove sequences less than 1000 nts.
(2) Sequences will be scanned for terminally redundant ends (circular sequences) and only sequences that are determined to be circular will be kept.
(3) Determines if circular contig has any ORFs of 80 AA or larger or else discards sequence.
(4)  (optional: see '-keep_known') Uses BLASTN against GenBank 'nt' to disregard any circular sequences that are >90% identical to known sequences across a >500bp window
(5) Uses Circlator ([](https://github.com/sanger-pathogens/circlator)) to rotate circular contigs so that a non-intragenic start codon of one of the ORFs will be the wrap point.

_-needs_rotation option_
(1) The specified .fasta file of pre-wrapped circular contigs will be used for Circlator (https://github.com/sanger-pathogens/circlator) to rotate circular contigs so that a non-intragenic start codon of one of the ORFs will be the wrap point.

_-given_linear or -given_circular options_
(1) given .fasta file will be kept more or less as is for annotation

### **Part 2: ORF Annotation and Taxonomy Guessing**
(1) Uses BLASTX against a custom virus + plasmid database (derived from GenBank 'nr') to guess family-level taxon of each circular sequence (note: this taxonomy method should not be considered definitive)
(2) Translates each ORF of 80 AA or larger
(3) Uses RPS-BLAST to predict function of each ORF by aligning to known NCBI Conserved Domains
(4) Generates a .tbl file of RPS-BLAST results
(5) Takes ORFs without RPS-BLAST hits and queries the GenBank viral database (predicted viruses) or nr database (predicted plasmids) with BLASTP
(6) Generates a .tbl file of BLASTP results
(7) Takes ORFs without any BLASTP hits and tries to find structural homology to known proteins using HHblits (databases: uniprot20, pdb70, scop70, pfam_31, NCBI_CD)
(8) Generates a .tbl file of HHblits results
(9) Removes ORFs within ORFs that do not contain conserved sequences
(10) Combines all .tbl files into a master .tbl file
(11) Generates a name for each virus/plasmid based on taxonomic results and nature of sample
(12) Generates properly formatted .fsa and .tbl files in a separate directory
(13) Uses tbl2asn to make gbf, val, and sqn files  

 Review all files and make any necessary changes before submitting to GenBank ([](https://www.ncbi.nlm.nih.gov/books/NBK53709/)). You will submit the .sqn file. Changes made to .gbf files will NOT be reflected on the .sqn file, unfortunately. However, the .gbf file should be looked at in a sequence map viewer such as MacVector, Geneious, or Ugene (free to use).

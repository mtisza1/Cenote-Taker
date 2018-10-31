#!/bin/bash

### This will only work in a Linux enviroment, I believe

# This is Cenote-Taker for Biowulf
# This version requires nucleotide fasta file of contigs (ideally from SPAdes with '--plasmid' option) as input. 


# Cenote-Taker Logo
echo "$(tput setaf 2)00000000000000000000000000"
echo "00000000000000000000000000"
echo "000000000$(tput setaf 4)^^^^^^^^$(tput setaf 2)000000000"
echo "000000$(tput setaf 4)^^^^^^^^^^^^^^$(tput setaf 2)000000"
echo "00000$(tput setaf 4)^^^^^$(tput setaf 3)CENOTE$(tput setaf 4)^^^^^$(tput setaf 2)00000"
echo "00000$(tput setaf 4)^^^^^$(tput setaf 3)TAKER!$(tput setaf 4)^^^^^$(tput setaf 2)00000"
echo "00000$(tput setaf 4)^^^^^^^^^^^^^^^^$(tput setaf 2)00000"
echo "000000$(tput setaf 4)^^^^^^^^^^^^^^$(tput setaf 2)000000"
echo "000000000$(tput setaf 4)^^^^^^^^$(tput setaf 2)000000000"
echo "00000000000000000000000000"
echo "00000000000000000000000000$(tput sgr 0)"

# Loading all the modules and environments for biowulf
### Regarding all the lines in this section: I don't know how your cloud server works but all these commands are specific to the NIH Biowulf system.
### 'source' and 'conda' commands are to open a python3 environment (as opposed to default python2 environment). Python3 is required for circlator.

source /data/tiszamj/conda/etc/profile.d/conda.sh
conda activate cenote_taker1
### find and install all of the following tools to your path.
### Here is the link to the 'edirect' tool: ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/

module load samtools || echo "$(tput setaf 4)unable to load samtools module $(tput sgr 0)"
module load mummer || echo "$(tput setaf 4)unable to load mummer module $(tput sgr 0)"
module load prodigal || echo "$(tput setaf 4)unable to load prodigal module $(tput sgr 0)" 
module load bioawk || echo "$(tput setaf 4)unable to load bioawk module $(tput sgr 0)"
module load blast || echo "$(tput setaf 4)unable to load blast module $(tput sgr 0)" 
module load hhsuite || echo "$(tput setaf 4)unable to load hhsuite module $(tput sgr 0)"
module load bwa || echo "$(tput setaf 4)unable to load BWA module $(tput sgr 0)"
module load edirect || echo "$(tput setaf 4)unable to load edirect module $(tput sgr 0)"
module load kronatools || echo "$(tput setaf 4)unable to load kronatools module $(tput sgr 0)"

# Setting input parameters
original_spades_contigs=$2
run_title=$3
isolation_source=$4
collection_date=$5
metagenome_type=$6
srr_number=$7
srx_number=$8
biosample=$9
bioproject=${10}
template_file=${11}


# Making output folder
if [ ! -d "$run_title" ]; then
	mkdir "$run_title"
fi 

# Taking contig input options
if  [[ $1 = "-given_circular" ]] ||  [[ $1 = "-given_linear" ]]; then
    if  [[ $1 = "-given_circular" ]]; then
    	echo "Option -given_circular turned on"
    elif  [[ $1 = "-given_linear" ]]; then
    	  echo "Option -given_linear turned on"
    fi
# splitting and renaming given contigs as fasta header ("ct_" + last 16 characters)
    if [ ${original_spades_contigs: -6} == ".fasta" ]; then
		echo "$(tput setaf 5)File with .fasta extension detected.$(tput sgr 0)"
### install bioawk
		bioawk -v run_var="$run_title" -c fastx '{ print ">"$name" "run_var NR; print $seq }' $original_spades_contigs > ${original_spades_contigs%.fasta}.rename.fasta ;
		cd $run_title
### csplit may or may not be in your path already,
		csplit -z ../${original_spades_contigs%.fasta}.rename.fasta '/>/' '{*}' --prefix=temp. --suffix-format=%03d.fasta;
	else
		echo "$(tput setaf 4)File with .fasta extension not detected as first input. Exiting.$(tput sgr 0)" ;
		exit
	fi

	for z in temp.*.fasta ; do 
		new_title=$( grep ">" $z | sed 's/>//g' | cut -d " " -f1 ) ; 
		echo $new_title ; 
		if [ ${#new_title} -ge 17 ]; then
			mv $z ct_${new_title: -16}.fasta ; 
		else
			mv $z ct_${new_title}.fasta ;
		fi
	done

	# Changing to output directory
	cd $run_title
	novel_fastas=$( ls *.fasta )

	if [ -z "$novel_fastas" ] ; then
		echo "$(tput setaf 4)No sequences detected in this directory. Exiting. $(tput sgr 0)" 
		exit
	else
		echo "$(tput setaf 5)Sequence(s) detected in directory$(tput sgr 0)"
		echo " "
	fi

	echo "$(tput setaf 5)Automatically annotating novel sequences:$(tput sgr 0)"
	for lineg in $novel_fastas ; do
		echo "$(tput setaf 5)"$lineg"$(tput sgr 0)"
	done
	echo " "

# script for -needs_rotation option
elif [[ $1 = "-needs_rotation" ]]; then
    echo "Option -needs_rotation turned on" 

# splitting and renaming given contigs as fasta header ("ct_" + last 16 characters)
    if [ ${original_spades_contigs: -6} == ".fasta" ]; then
		echo "$(tput setaf 5)File with .fasta extension detected.$(tput sgr 0)"
		bioawk -v run_var="$run_title" -c fastx '{ print ">"$name" "run_var NR; print $seq }' $original_spades_contigs > ${original_spades_contigs%.fasta}.rename.fasta ;
		cd $run_title
		csplit -z ../${original_spades_contigs%.fasta}.rename.fasta '/>/' '{*}' --prefix=temp. --suffix-format=%03d.fasta;
	else
		echo "$(tput setaf 4)File with .fasta extension not detected as first input. Exiting.$(tput sgr 0)" ;
		exit
	fi

	for z in temp.*.fasta ; do 
		new_title=$( grep ">" $z | sed 's/>//g' | cut -d " " -f1 ) ; 
		echo $new_title ; 
		if [ ${#new_title} -ge 17 ]; then
			mv $z ct_${new_title: -16}.fasta ; 
		else
			mv $z ct_${new_title}.fasta ;
		fi
	done

	# Changing to output directory
	cd $run_title
	novel_fastas=$( ls *.fasta )

	if [ -z "$novel_fastas" ] ; then
		echo "$(tput setaf 4)No sequences detected in this directory. Exiting. $(tput sgr 0)" 
		exit
	else
		echo "$(tput setaf 5)Sequence(s) detected in directory$(tput sgr 0)"
		echo " "
	fi

	echo "$(tput setaf 5)Automatically annotating novel sequences:$(tput sgr 0)"
	for lineg in $novel_fastas ; do
		echo "$(tput setaf 5)"$lineg"$(tput sgr 0)"
	done
	echo " "

	# Rotating each sequence to put a non-intragenic start codon as the first basepair of the contig
	for nucl_fa in $novel_fastas ; do
	echo "$(tput setaf 5)rotating "$nucl_fa" to put an ORF at beginning of sequence so that no ORFs overlap the breakpoint $(tput sgr 0)"
### for the the getorf function, install the EMBOSS suite of tools
	/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/getorf -circular -minsize 240 -find 3 -sequence $nucl_fa -outseq ${nucl_fa%.fasta}.nucl_orfs.fa ; 

	grep ">" ${nucl_fa%.fasta}.nucl_orfs.fa > ${nucl_fa%.fasta}.nucl_orfs.txt
	cat "${nucl_fa%.fasta}.nucl_orfs.txt" | while read liner ; do
	start_base=$( echo $liner | sed 's/.*\[\(.*\) -.*/\1/' )
	end_base=$( echo $liner | sed 's/.*- \(.*\)\].*/\1/' )
	length=$(( $start_base-$end_base ))
	abso_length=$( echo $length | sed 's/-//g' )
	if [ $abso_length -gt 239 ]; then
		if [[ "$end_base" -gt "$start_base" ]]; then
			for ((counter_f=(( $start_base + 1 ));counter_f<=(( $end_base + 3 ));counter_f++)); do
				echo "$counter_f" >> ${nucl_fa%.fasta}.bad_starts.txt
				
			done
		elif [[ "$start_base" -gt "$end_base" ]]; then
			for ((counter_r=(( $end_base - 3 ));counter_r<=(( $start_base - 1 ));counter_r++)) ; do
				echo "$counter_r" >> ${nucl_fa%.fasta}.bad_starts.txt
			done
		fi
	fi

	done

	cat "${nucl_fa%.fasta}.nucl_orfs.txt" | while read liner ; do
		starter_base=$( echo $liner | sed 's/.*\[\(.*\) -.*/\1/' )
		if grep -q "$starter_base" ${nucl_fa%.fasta}.bad_starts.txt ; then
			continue
		else
			echo $liner >> ${nucl_fa%.fasta}.good_start_orfs.txt
		fi

	done
	if [ -s "${nucl_fa%.fasta}.good_start_orfs.txt" ]; then	
		cut -d " " -f1 ${nucl_fa%.fasta}.good_start_orfs.txt | head -n1 | sed 's/>//g' > ${nucl_fa%.fasta}.starting_orf.txt
		bioawk -c fastx '{ print $name, $seq, length($seq) }' ${nucl_fa%.fasta}.nucl_orfs.fa | grep -f ${nucl_fa%.fasta}.starting_orf.txt | sed '/--/d' | head -n1 | awk '{print ">"$1, $3; print $2}' > ${nucl_fa%.fasta}.starting_orf.1.fa ;
### After installing all of its dependencies all of which are mentioned at the top of this script, install circlator
		circlator fixstart --genes_fa ${nucl_fa%.fasta}.starting_orf.1.fa $nucl_fa ${nucl_fa%.fasta}.rotate ;
	else
		cp $nucl_fa ${nucl_fa%.fasta}.no_100AA_ORFs.fasta
	fi
	done 

# script for -default option
elif [[ $1 = "-default" ]]; then
    echo "Option -default turned on" 
    # Removing contigs under 1000bp and detecting circular contigs
	if [ ${original_spades_contigs: -6} == ".fasta" ]; then
		echo "$(tput setaf 5)File with .fasta extension detected, attempting to keep contigs over 1000bp and find circular sequences with apc.pl$(tput sgr 0)"
		bioawk -v run_var="$run_title" -c fastx '{ if(length($seq) > 1000) { print ">"run_var NR" "$name; print $seq }}' $original_spades_contigs > ${original_spades_contigs%.fasta}.over_1000.fasta ;
		cd $run_title
### You'll need 'last' from: http://last.cbrc.jp/ 
### go into the apc_biowulf1.pl script and change the path to lastal and lastdb
		perl /data/tiszamj/mike_tisza/auto_annotation_pipeline/apc_biowulf1.pl -b $run_title ../${original_spades_contigs%.fasta}.over_1000.fasta ;
		rm apc_aln*
		for permu_file in permuted.*.fa ; do 
			mv $permu_file "$run_title"${permu_file#permuted.} ; 
			for fa1 in $run_title*.fa ; do 
				mv $fa1 $run_title${fa1#$run_title.}sta ; 
			done 
		done
	else
		echo "$(tput setaf 4)File with .fasta extension not detected as first input. Exiting.$(tput sgr 0)" ;
		exit
	fi

	# Changing to output directory
	cd $run_title


	# Detecting whether any circular contigs were present
	original_fastas=$( ls *.fasta )
	# "$(tput setaf 5)$var1$(tput sgr 0)"
	if [ -z "$original_fastas" ] ; then
		echo "$(tput setaf 4)No circular fasta files detected. Exiting. $(tput sgr 0)" 
		exit
	else
		echo "$(tput setaf 5)Circular fasta file(s) detected$(tput sgr 0)"
		echo " "
	fi

	# Checking whether any circular contigs are >90% identical to any sequence in NCBI nt database using BLASTN. If so, disregarding.
	if  [[ ${12} = "-keep_known" ]]; then
		echo " "
	else
		echo "$(tput setaf 5)Checking Novelty of:$(tput sgr 0)"
		for line in $original_fastas ; do
			echo "$(tput setaf 5)"$line"$(tput sgr 0)"
		done
	fi
	echo " "

	for circle in $original_fastas ; do
	if  [[ ${12} = "-keep_known" ]]; then
		echo " -keep_known option used. Not searching for or disregarding known sequences"
	else
### install the blast suite and download and format the 'nt' database from genbank.
		blastn -db /fdb/blastdb/nt -query $circle -evalue 1e-50 -num_threads 56 -outfmt "6 qseqid sseqid stitle pident length" -qcov_hsp_perc 50 -num_alignments 3 -out ${circle%.fasta}.blastn.out ;
		if [ -s "${circle%.fasta}.blastn.out" ]; then
			sed 's/ /-/g' ${circle%.fasta}.blastn.out | awk '{if ($4 > 90) print}' | awk '{if ($5 > 500) print}' > ${circle%.fasta}.blastn.notnew.out ;
		fi
		if [ -s "${circle%.fasta}.blastn.notnew.out" ]; then
			echo "$(tput setaf 4)"$circle" is not a novel species (>90% identical to sequence in nt database) and will not be annotated.$(tput sgr 0)"
			mv $circle ${circle%.fasta}.notnew.fna ;
		else
			echo "$(tput setaf 5)"$circle" appears to be a novel sequence (no close (>90% nucleotide) matches to sequences in nt database).$(tput sgr 0)"
		fi
	fi
	done

	rm *.blastn.out

	# echoing novel circular sequences
	novel_fastas=$( ls *.fasta )
	if  [[ ${12} = "-keep_known" ]]; then
		echo "$(tput setaf 5)Automatically annotating sequences:$(tput sgr 0)"
		for lineg in $novel_fastas ; do
			echo "$(tput setaf 5)"$lineg"$(tput sgr 0)"
		done
	else
		echo " "
		if [ -z "$novel_fastas" ] ; then
			echo "$(tput setaf 4)No novel circular sequences detected. Exiting. $(tput sgr 0)" 
			exit
		else
			echo "$(tput setaf 5)Novel circular sequence(s) detected$(tput sgr 0)"
			echo " "
		fi

		echo "$(tput setaf 5)Automatically annotating novel sequences:$(tput sgr 0)"
		for lineg in $novel_fastas ; do
			echo "$(tput setaf 5)"$lineg"$(tput sgr 0)"
		done
		echo " "
	fi

	# Rotating each sequence to put a non-intragenic start codon as the first basepair of the contig
	for nucl_fa in $novel_fastas ; do
	echo "$(tput setaf 5)rotating "$nucl_fa" to put an ORF at beginning of sequence so that no ORFs overlap the breakpoint $(tput sgr 0)"
	/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/getorf -circular -minsize 240 -find 3 -sequence $nucl_fa -outseq ${nucl_fa%.fasta}.nucl_orfs.fa ; 

	grep ">" ${nucl_fa%.fasta}.nucl_orfs.fa > ${nucl_fa%.fasta}.nucl_orfs.txt
	cat "${nucl_fa%.fasta}.nucl_orfs.txt" | while read liner ; do
	start_base=$( echo $liner | sed 's/.*\[\(.*\) -.*/\1/' )
	end_base=$( echo $liner | sed 's/.*- \(.*\)\].*/\1/' )
	length=$(( $start_base-$end_base ))
	abso_length=$( echo $length | sed 's/-//g' )
	if [ $abso_length -gt 239 ]; then
		if [[ "$end_base" -gt "$start_base" ]]; then
			for ((counter_f=(( $start_base + 1 ));counter_f<=(( $end_base + 3 ));counter_f++)); do
				echo "$counter_f" >> ${nucl_fa%.fasta}.bad_starts.txt
				
			done
		elif [[ "$start_base" -gt "$end_base" ]]; then
			for ((counter_r=(( $end_base - 3 ));counter_r<=(( $start_base - 1 ));counter_r++)) ; do
				echo "$counter_r" >> ${nucl_fa%.fasta}.bad_starts.txt
			done
		fi
	fi

	done

	cat "${nucl_fa%.fasta}.nucl_orfs.txt" | while read liner ; do
		starter_base=$( echo $liner | sed 's/.*\[\(.*\) -.*/\1/' )
		if grep -q "$starter_base" ${nucl_fa%.fasta}.bad_starts.txt ; then
			continue
		else
			echo $liner >> ${nucl_fa%.fasta}.good_start_orfs.txt
		fi

	done
	if [ -s "${nucl_fa%.fasta}.good_start_orfs.txt" ]; then	
		cut -d " " -f1 ${nucl_fa%.fasta}.good_start_orfs.txt | head -n1 | sed 's/>//g' > ${nucl_fa%.fasta}.starting_orf.txt
		bioawk -c fastx '{ print $name, $seq, length($seq) }' ${nucl_fa%.fasta}.nucl_orfs.fa | grep -f ${nucl_fa%.fasta}.starting_orf.txt | sed '/--/d' | head -n1 | awk '{print ">"$1, $3; print $2}' > ${nucl_fa%.fasta}.starting_orf.1.fa ;
		circlator fixstart --genes_fa ${nucl_fa%.fasta}.starting_orf.1.fa $nucl_fa ${nucl_fa%.fasta}.rotate ;
	else
		cp $nucl_fa ${nucl_fa%.fasta}.no_100AA_ORFs.fasta
	fi
	done 
else
    echo "$(tput setaf 5)You did not use any contig input option. Please rerun with option.$(tput sgr 0)"
    exit
fi

############
if  [[ $1 = "-given_circular" ]] ||  [[ $1 = "-given_linear" ]]; then
	for nucl_fa in $novel_fastas ; do
		cp $nucl_fa ${nucl_fa%.fasta}.rotate.fasta
	done
fi

# Performing BLASTX of each contig against database of viral and plasmid proteins to guess taxonomy
for nucl_fa in $novel_fastas ; do
if [ -s "${nucl_fa%.fasta}.rotate.fasta" ]; then
	echo "$(tput setaf 5)Guessing taxonomy for sequence "${nucl_fa%.fasta}.rotate.fasta" by BLASTX against virus and plasmid protein database.$(tput sgr 0)"
### I have a custom-ish database of viral and plasmid proteins. The viral proteins are directly taken from refseq. The plasmid proteins were pulled out of the 'nr' bacterial database by me. I will send this database. 
### You may choose to use a slightly different database, but for this step, don't use a database of viruses you've annotated because downstream parts of cenote-taker need the genbank taxonomy info.
	blastx -evalue 1e-4 -outfmt "6 qseqid stitle pident evalue length" -num_threads 56 -num_alignments 1 -db /data/tiszamj/mike_tisza/auto_annotation_pipeline/blast_DBs/virus_and_plasmid_proteins -query ${nucl_fa%.fasta}.rotate.fasta -out ${nucl_fa%.fasta}.tax_guide.blastx.out ;
	if [ ! -s "${nucl_fa%.fasta}.tax_guide.blastx.out" ]; then
		echo "No homologues found" > ${nucl_fa%.fasta}.tax_guide.blastx.out ;
	else
		echo "$(tput setaf 5)"$nucl_fa" likely represents a novel virus or plasmid. Getting hierarchical taxonomy info.$(tput sgr 0)"
### ktClassifyBLAST is available when you install kronatools
		ktClassifyBLAST -o ${nucl_fa%.fasta}.tax_guide.blastx.tab ${nucl_fa%.fasta}.tax_guide.blastx.out
		taxid=$( tail -n1 ${nucl_fa%.fasta}.tax_guide.blastx.tab | cut -f2 )
### efetch is from edirect, as is xtract, I believe.
		efetch -db taxonomy -id $taxid -format xml | /data/tiszamj/mike_tisza/xtract.Linux -pattern Taxon -element Lineage >> ${nucl_fa%.fasta}.tax_guide.blastx.out
	fi
else
	echo "$(tput setaf 4)"$nucl_fa" could not be rotated. Likely there were no ORFs of at least 100AA.$(tput sgr 0)" 

fi
done


# Extracting ORFs >240bp, (>180bp for inoviruses/plasmids)
for nucl_fa in $novel_fastas ; do
if [ -s "${nucl_fa%.fasta}.rotate.fasta" ]; then
	echo "$(tput setaf 5)"$nucl_fa" taxonomy guessed. Continuing to ORF translation...$(tput sgr 0)"

	if [[ $1 = "-given_linear" ]]; then
		if grep -q "Inovir" ${nucl_fa%.fasta}.tax_guide.blastx.out || grep -q "plasmid" ${nucl_fa%.fasta}.tax_guide.blastx.out ; then
			/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/getorf -find 1 -minsize 180 -sequence ${nucl_fa%.fasta}.rotate.fasta -outseq ${nucl_fa%.fasta}.rotate.AA.fasta ;
		else
			/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/getorf -find 1 -minsize 240 -sequence ${nucl_fa%.fasta}.rotate.fasta -outseq ${nucl_fa%.fasta}.rotate.AA.fasta ;
		fi
	else
		if grep -q "Inovir" ${nucl_fa%.fasta}.tax_guide.blastx.out ; then
			/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/getorf -find 1 -minsize 180 -sequence ${nucl_fa%.fasta}.rotate.fasta -outseq ${nucl_fa%.fasta}.rotate.AA.fasta ;
		else
			/data/tiszamj/mike_tisza/EMBOSS-6.6.0/emboss/getorf -find 1 -minsize 240 -sequence ${nucl_fa%.fasta}.rotate.fasta -outseq ${nucl_fa%.fasta}.rotate.AA.fasta ;
		fi
	fi


	bioawk -c fastx '{FS="\t"; OFS=" "} {print ">"$name $3, $4, $5, $6, $7; print $seq}' ${nucl_fa%.fasta}.rotate.AA.fasta > ${nucl_fa%.fasta}.rotate.AA.sorted.fasta ;
fi
done

# Conducting RPS-BLAST against CDD on translated ORFs
for nucl_fa in $novel_fastas ; do
if [ -s "${nucl_fa%.fasta}.rotate.AA.sorted.fasta" ]; then

	echo "$(tput setaf 5)"$nucl_fa" Continuing to RPS-BLAST NCBI CDD domains database for each ORF...$(tput sgr 0)" 
### CDD database can be found here: ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd
	rpsblast -evalue 1e-4 -num_descriptions 5 -num_threads 56 -line_length 100 -num_alignments 1 -db /data/tiszamj/mike_tisza/auto_annotation_pipeline/cdd_rps_db/Cdd -query ${nucl_fa%.fasta}.rotate.AA.sorted.fasta -out ${nucl_fa%.fasta}.rotate.AA.rpsblast.out ;
	echo "$(tput setaf 5)RPS-BLAST of "${nucl_fa%.fasta}.rotate.AA.sorted.fasta" complete.$(tput sgr 0)"
	echo " "
else
	echo "$(tput setaf 4)ORF file for "$nucl_fa" is empty. This circle might have no ORFS over 100AA.$(tput sgr 0)"
	echo " "
fi
done
#rm ${nucl_fa%.fasta}.nucl_orfs.fa ${nucl_fa%.fasta}.rotate.detailed.log ${nucl_fa%.fasta}.rotate.log ${nucl_fa%.fasta}.rotate.promer.promer ${nucl_fa%.fasta}.rotate.promer.contigs_with_ends.fa ${nucl_fa%.fasta}.rotate.prodigal.for_prodigal.fa ${nucl_fa%.fasta}.rotate.prodigal.prodigal.gff;


# Generating tbl file from RPS-BLAST results
### set path to this perl script. 
perl /data/tiszamj/mike_tisza/auto_annotation_pipeline/rpsblastreport2tbl_mt_annotation_pipe_biowulf.pl ;

for nucl_fa in $novel_fastas ; do 
if [ -s "${nucl_fa%.fasta}.NT.tbl" ]; then
	echo "$(tput setaf 5)"$nucl_fa" tbl made from RPS-BLAST hits...$(tput sgr 0)"
else
	echo "$(tput setaf 4) RPS-BLAST tbl for "$nucl_fa" not detected.$(tput sgr 0)"
fi
done

# Conducting BLASTP on ORFs unrecognized by RPS-BLAST (nr virus or nr all database)
for feat_tbl1 in *.NT.tbl ; do
grep -e 'hypothetical protein' -e 'unnamed protein product' -e 'predicted protein' -e 'Uncharacterized protein' -e 'Domain of unknown function' -B2 $feat_tbl1 | grep "^[0-9]" | awk '{print $1 " - " $2}' > ${feat_tbl1%.NT.tbl}.for_hhpred.txt ;
grep -f ${feat_tbl1%.NT.tbl}.for_hhpred.txt -A1 ${feat_tbl1%.NT.tbl}.rotate.AA.sorted.fasta | sed '/--/d' > ${feat_tbl1%.NT.tbl}.rotate.rps_nohits.fasta ;
if grep -q "(plasmid)" ${feat_tbl1%.NT.tbl}.tax_guide.blastx.out ; then
	if [ -s "${feat_tbl1%.NT.tbl}.rotate.rps_nohits.fasta" ]; then
		echo "$(tput setaf 5)"$nucl_fa" is likely a plasmid... Continuing to BLASTP NCBI nr database for each ORF that had no hits in CDD...$(tput sgr 0)" 
### I use both the 'nr' and 'nr/viral' databases separately. Available from genbank.
		blastp -evalue 1e-4 -num_descriptions 5 -num_threads 56 -num_alignments 1 -db /fdb/blastdb/nr -query ${feat_tbl1%.NT.tbl}.rotate.rps_nohits.fasta -out ${feat_tbl1%.NT.tbl}.rotate.blastp.out ;
		echo "$(tput setaf 5)BLASTP of "${feat_tbl1%.NT.tbl}.rotate.rps_nohits.fasta" complete.$(tput sgr 0)"
	fi
else
	if [ -s "${feat_tbl1%.NT.tbl}.rotate.rps_nohits.fasta" ]; then
		echo "$(tput setaf 5)"$nucl_fa" is likely a virus... Continuing to BLASTP NCBI virus database for each ORF that had no hits in CDD...$(tput sgr 0)" 
		blastp -evalue 1e-4 -num_descriptions 5 -num_threads 56 -num_alignments 1 -db /fdb/blastdb/viral -query ${feat_tbl1%.NT.tbl}.rotate.rps_nohits.fasta -out ${feat_tbl1%.NT.tbl}.rotate.blastp.out ;
		echo "$(tput setaf 5)BLASTP of "${feat_tbl1%.NT.tbl}.rotate.rps_nohits.fasta" complete.$(tput sgr 0)"
	fi
fi
done

# Generating tbl file from BLASTP results
### set path to this perl script. 
perl /data/tiszamj/mike_tisza/auto_annotation_pipeline/blastpreport2tbl_mt_annotation_pipe_biowulf2.pl ;

for feat_tbl1 in *.NT.tbl ; do
if [ -s "${feat_tbl1%.NT.tbl}.BLASTP.tbl" ]; then
	echo "$(tput setaf 5)"${feat_tbl1%.NT.tbl}": tbl made from BLASTP hits. Splitting fasta files for HHBLITS...$(tput sgr 0)"
else
	echo "$(tput setaf 4) BLASTP tbl for "${feat_tbl1%.NT.tbl}" not detected.$(tput sgr 0)"
fi
done

# Grabbing ORFs wihout BLASTP hits and separating them into individual files for HHBlits
for blastp_tbl1 in *.BLASTP.tbl ; do
grep -e 'hypothetical protein' -e 'unnamed protein product' -e 'predicted protein' -e 'Uncharacterized protein' -e 'Uncharacterized conserved protein' -e 'unknown' -e 'Uncharacterised protein' -B2 $blastp_tbl1 | grep "^[0-9]" | awk '{print $1 " - " $2}' > ${blastp_tbl1%.BLASTP.tbl}.for_hhpred.txt ;
grep -f ${blastp_tbl1%.BLASTP.tbl}.for_hhpred.txt -A1 ${blastp_tbl1%.BLASTP.tbl}.rotate.AA.sorted.fasta | sed '/--/d' > ${blastp_tbl1%.BLASTP.tbl}.rotate.blast_hypo.fasta ;
csplit -z ${blastp_tbl1%.BLASTP.tbl}.rotate.blast_hypo.fasta '/>/' '{*}' --prefix=${blastp_tbl1%.BLASTP.tbl}.rotate --suffix-format=%02d.for_hhpred.fasta; 
done

# Running HHBlits on remaining ORFs
dark_orf_list=$( ls *.for_hhpred.fasta )

for dark_orf in $dark_orf_list ; do
echo "$(tput setaf 5)Running HHBLITS on "$dark_orf" now.$(tput sgr 0)"
### HHBLITS should be available here: https://github.com/soedinglab/hh-suite
### most of the databases are here: http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/
### However, the CDD database was sent to me specially by the people at Tuebingen.
hhblits -i $dark_orf -d /fdb/hhsuite/uniprot20_2016_02/uniprot20_2016_02 -d /data/tiszamj/mike_tisza/auto_annotation_pipeline/pdb70/pdb70 -d /data/tiszamj/mike_tisza/auto_annotation_pipeline/scop70/scop70_1.75 -d /data/tiszamj/mike_tisza/auto_annotation_pipeline/pfam_31_db/pfam -d /data/tiszamj/mike_tisza/auto_annotation_pipeline/cdd_db/NCBI_CD/NCBI_CD -o ${dark_orf%.for_hhpred.fasta}.out.hhr -cov 50 -cpu 56 -maxmem 70 -p 90 -Z 20 -z 0 -b 0 -B 10 -ssm 2 -sc 1 -E 1  ;
cat ${dark_orf%.for_hhpred.fasta}.out.hhr >> ${dark_orf%.rotate*.for_hhpred.fasta}.rotate.out_all.hhr ;
done

rm *.rotate.AA.fasta

# Generating tbl file from HHBlits results
### set path to this perl script. 
perl /data/tiszamj/mike_tisza/auto_annotation_pipeline/hhpredreport2tbl_mt_annotation_pipe_biowulf1_gjs_edits.pl ;

for HH_tbl1 in *.HH.tbl ; do 
sed 's/OS=.*//g; s/ ;//g; s/UniProtKB:>\([0-9][A-Z].*\)/PDB:\1/g; s/UniProtKB:>tr|.*|\(.\)/UniProtKB:\1/g; s/UniProtKB:>\([a-z].*\)/Scop:\1/g; s/UniProtKB:>\(PF.*\)/PFAM:\1/g; s/ is .*//g; s/ are .*//g' $HH_tbl1 > ${HH_tbl1%.HH.tbl}.HH2.tbl
done

# Combining tbl files from all search results

for feat_tbl2 in *.NT.tbl ; do 
if [ -s "${feat_tbl2%.NT.tbl}.HH2.tbl" ]; then
	grep -v -e 'hypothetical protein' -e 'unnamed protein product' -e 'Predicted protein' -e 'predicted protein' -e 'Uncharacterized protein' -e 'Uncharacterized conserved protein' -e 'unknown' -e 'Uncharacterised protein' ${feat_tbl2%.NT.tbl}.BLASTP.tbl | grep -A1 -B2 'product' | grep -v ">Feature" | sed '/--/d' > ${feat_tbl2%.NT.tbl}.tmp.tbl ; 
	grep -v -e 'hypothetical protein' -e 'unnamed protein product' -e 'Predicted protein' -e 'predicted protein' -e 'Uncharacterized protein' -e 'Domain of unknown function' $feat_tbl2 | grep -A1 -B2 'product' | grep -v ">Feature" | sed '/--/d' >> ${feat_tbl2%.NT.tbl}.tmp.tbl ; 
	cat ${feat_tbl2%.NT.tbl}.HH2.tbl <(echo) ${feat_tbl2%.NT.tbl}.tmp.tbl > ${feat_tbl2%.NT.tbl}.comb.tbl ; 
elif [ -s "${feat_tbl2%.NT.tbl}.BLASTP.tbl" ]; then
	grep -v -e 'hypothetical protein' -e 'unnamed protein product' -e 'Predicted protein' -e 'predicted protein' -e 'Uncharacterized protein' -e 'Domain of unknown function' $feat_tbl2 | grep -A1 -B2 'product' | grep -v ">Feature" | sed '/--/d' > ${feat_tbl2%.NT.tbl}.tmp.tbl ; 
	cat ${feat_tbl2%.NT.tbl}.BLASTP.tbl <(echo) ${feat_tbl2%.NT.tbl}.tmp.tbl > ${feat_tbl2%.NT.tbl}.BLASTPcomb.tbl ; 
fi
done

rm *tmp.tbl

#removing hypothetical ORFs-within-ORFs
for feat_tbl3 in *.NT.tbl ; do
	grep "^[0-9]" $feat_tbl3 | awk '{FS="\t"; OFS="\t"} {print $1, $2}' > ${feat_tbl3%.NT.tbl}.all_start_stop.txt ;
	cat "${feat_tbl3%.NT.tbl}.all_start_stop.txt" | while read linev ; do
		all_start=$( echo $linev | cut -d " " -f1 )
		all_end=$( echo $linev | cut -d " " -f2 )
		if [[ "$all_end" -gt "$all_start" ]]; then
			for ((counter_f=(( $all_start + 1 ));counter_f<=$all_end;counter_f++)); do
				echo " " "$counter_f" " " >> ${feat_tbl3%.NT.tbl}.used_positions.txt
				
			done
		elif [[ "$all_start" -gt "$all_end" ]]; then
			for ((counter_r=$all_end;counter_r<=(( $all_start - 1 ));counter_r++)) ; do
				echo " " "$counter_r" " " >> ${feat_tbl3%.NT.tbl}.used_positions.txt
			done
		fi
	done

if [ -s "${feat_tbl3%.NT.tbl}.comb.tbl" ]; then
	sed 's/(Fragment)//g; s/\. .*//g; s/{.*//g; s/\[.*//g; s/Putative hypothetical protein/hypothetical protein/g; s/Uncultured bacteri.*/hypothetical protein/g; s/RNA helicase$/helicase/g; s/Os.*/hypothetical protein/g; s/\.$//g; s/Unplaced genomic scaffold.*/hypothetical protein/g; s/Putative hypothetical protein/hypothetical protein/g; s/Contig.*/hypothetical protein/g; s/Uncharacterized protein/hypothetical protein/g; s/uncharacterized protein/hypothetical protein/g; s/Uncharacterised protein/hypothetical protein/g' ${feat_tbl3%.NT.tbl}.comb.tbl | sed '/--/d' > ${feat_tbl3%.NT.tbl}.comb2.tbl ; 
	grep -e 'hypothetical protein' -B2 ${feat_tbl3%.NT.tbl}.comb2.tbl | grep "^[0-9]" | awk '{FS="\t"; OFS="\t"} {print $1, $2}' > ${feat_tbl3%.NT.tbl}.hypo_start_stop.txt ;

	cat "${feat_tbl3%.NT.tbl}.hypo_start_stop.txt" | while read liney ; do
		loc_start=$( echo $liney | cut -d " " -f1 )
		loc_end=$( echo $liney | cut -d " " -f2 )
		loc1_start=$( echo " " "$loc_start" " ")
	if grep -q "$loc1_start" ${feat_tbl3%.NT.tbl}.used_positions.txt ; then 
		echo "$loc1_start"
		if [[ "$loc_end" -gt "$loc_start" ]]; then
			f_end=$(( $loc_end + 1 ))
			f1_end=$( echo " " "$f_end" " ")
			echo "$f1_end"
			if grep -q "$f1_end" ${feat_tbl3%.NT.tbl}.used_positions.txt ; then
				echo "$loc_end" "end"
				echo "$liney" >> ${feat_tbl3%.NT.tbl}.remove_hypo.txt
			fi
		else
			r_end=$(( $loc_end - 1 ))
			r1_end=$( echo " " "$r_end" " ")

			echo "$r1_end"
			if grep -q "$r1_end" ${feat_tbl3%.NT.tbl}.used_positions.txt ; then
				echo "$loc_end" "end"
				echo "$liney" >> ${feat_tbl3%.NT.tbl}.remove_hypo.txt
			fi
		fi
	fi
	done
	if [ -s "${feat_tbl3%.NT.tbl}.remove_hypo.txt" ]; then
		grep ">Feature" ${feat_tbl3%.NT.tbl}.comb2.tbl | sed '/--/d' > ${feat_tbl3%.NT.tbl}.HH_BLAST.tbl
		grep -v -f ${feat_tbl3%.NT.tbl}.remove_hypo.txt ${feat_tbl3%.NT.tbl}.comb2.tbl | grep "^[0-9]" -A3 | sed '/--/d' >> ${feat_tbl3%.NT.tbl}.HH_BLAST.tbl ;
	else
		cp ${feat_tbl3%.NT.tbl}.comb2.tbl ${feat_tbl3%.NT.tbl}.HH_BLAST.tbl
	fi
elif [ -s "${feat_tbl3%.NT.tbl}.BLASTPcomb.tbl" ]; then
	sed 's/(Fragment)//g; s/\. .*//g; s/{.*//g; s/\[.*//g; s/Putative hypothetical protein/hypothetical protein/g; s/Uncultured bacteri.*/hypothetical protein/g; s/RNA helicase$/helicase/g; s/Os.*/hypothetical protein/g; s/\.$//g; s/Unplaced genomic scaffold.*/hypothetical protein/g; s/Putative hypothetical protein/hypothetical protein/g; s/Contig.*/hypothetical protein/g; s/Uncharacterized protein/hypothetical protein/g; s/uncharacterized protein/hypothetical protein/g; s/Uncharacterised protein/hypothetical protein/g' ${feat_tbl3%.NT.tbl}.BLASTPcomb.tbl | sed '/--/d' > ${feat_tbl3%.NT.tbl}.BLASTP2.tbl ; 
	grep -e 'hypothetical protein' -B2 ${feat_tbl3%.NT.tbl}.BLASTP2.tbl | grep "^[0-9]" | awk '{FS="\t"; OFS="\t"} {print $1, $2}' > ${feat_tbl3%.NT.tbl}.hypo_start_stop.txt ;
	cat "${feat_tbl3%.NT.tbl}.hypo_start_stop.txt" | while read liney ; do
		loc_start=$( echo $liney | cut -d " " -f1 )
		loc_end=$( echo $liney | cut -d " " -f2 )
		loc1_start=$( echo " " "$loc_start" " ")
	if grep -q "$loc1_start" ${feat_tbl3%.NT.tbl}.used_positions.txt ; then 
		echo "$loc1_start"
		if [[ "$loc_end" -gt "$loc_start" ]]; then
			f_end=$(( $loc_end + 1 ))
			f1_end=$( echo " " "$f_end" " ")
			echo "$f1_end"
			if grep -q "$f1_end" ${feat_tbl3%.NT.tbl}.used_positions.txt ; then
				echo "$loc_end" "end"
				echo "$liney" >> ${feat_tbl3%.NT.tbl}.remove_hypo.txt
			fi
		else
			r_end=$(( $loc_end - 1 ))
			r1_end=$( echo " " "$r_end" " ")

			echo "$r1_end"
			if grep -q "$r1_end" ${feat_tbl3%.NT.tbl}.used_positions.txt ; then
				echo "$loc_end" "end"
				echo "$liney" >> ${feat_tbl3%.NT.tbl}.remove_hypo.txt
			fi
		fi
	fi
	done
	if [ -s "${feat_tbl3%.NT.tbl}.remove_hypo.txt" ]; then
		grep ">Feature" ${feat_tbl3%.NT.tbl}.BLASTP2.tbl | sed '/--/d' > ${feat_tbl3%.NT.tbl}.BLASTP_final.tbl
		grep -v -f ${feat_tbl3%.NT.tbl}.remove_hypo.txt ${feat_tbl3%.NT.tbl}.BLASTP2.tbl | grep "^[0-9]" -A3 | sed '/--/d' >> ${feat_tbl3%.NT.tbl}.BLASTP_final.tbl
	else
		cp ${feat_tbl3%.NT.tbl}.BLASTP2.tbl ${feat_tbl3%.NT.tbl}.BLASTP_final.tbl
	fi
else
	sed 's/(Fragment)//g; s/\. .*//g; s/{.*//g; s/\[.*//g; s/Putative hypothetical protein/hypothetical protein/g; s/Uncultured bacteri.*/hypothetical protein/g; s/RNA helicase$/helicase/g; s/Os.*/hypothetical protein/g; s/\.$//g; s/Unplaced genomic scaffold.*/hypothetical protein/g; s/Putative hypothetical protein/hypothetical protein/g; s/Contig.*/hypothetical protein/g; s/Uncharacterized protein/hypothetical protein/g; s/uncharacterized protein/hypothetical protein/g; s/Uncharacterised protein/hypothetical protein/g' $feat_tbl3 | sed '/--/d' > ${feat_tbl3%.NT.tbl}.NT2.tbl ; 
	grep -e 'hypothetical protein' -B2 ${feat_tbl3%.NT.tbl}.NT2.tbl | grep "^[0-9]" | awk '{FS="\t"; OFS="\t"} {print $1, $2}' > ${feat_tbl3%.NT.tbl}.hypo_start_stop.txt ;
	cat "${feat_tbl3%.NT.tbl}.hypo_start_stop.txt" | while read liney ; do
		loc_start=$( echo $liney | cut -d " " -f1 )
		loc_end=$( echo $liney | cut -d " " -f2 )
		loc1_start=$( echo " " "$loc_start" " ")
	if grep -q "$loc1_start" ${feat_tbl3%.NT.tbl}.used_positions.txt ; then 
		echo "$loc1_start"
		if [[ "$loc_end" -gt "$loc_start" ]]; then
			f_end=$(( $loc_end + 1 ))
			f1_end=$( echo " " "$f_end" " ")
			echo "$f1_end"
			if grep -q "$f1_end" ${feat_tbl3%.NT.tbl}.used_positions.txt ; then
				echo "$loc_end" "end"
				echo "$liney" >> ${feat_tbl3%.NT.tbl}.remove_hypo.txt
			fi
		else
			r_end=$(( $loc_end - 1 ))
			r1_end=$( echo " " "$r_end" " ")

			echo "$r1_end"
			if grep -q "$r1_end" ${feat_tbl3%.NT.tbl}.used_positions.txt ; then
				echo "$loc_end" "end"
				echo "$liney" >> ${feat_tbl3%.NT.tbl}.remove_hypo.txt
			fi
		fi
	fi
	done
	if [ -s "${feat_tbl3%.NT.tbl}.remove_hypo.txt" ]; then
		grep ">Feature" ${feat_tbl3%.NT.tbl}.NT2.tbl | sed '/--/d' > ${feat_tbl3%.NT.tbl}.NT_final.tbl
		grep -v -f ${feat_tbl3%.NT.tbl}.remove_hypo.txt ${feat_tbl3%.NT.tbl}.NT2.tbl | grep "^[0-9]" -A3 | sed '/--/d' >> ${feat_tbl3%.NT.tbl}.NT_final.tbl
	else
		cp ${feat_tbl3%.NT.tbl}.NT2.tbl ${feat_tbl3%.NT.tbl}.NT_final.tbl
	fi
fi
done



# Making directory for sequin generation
if [ ! -d "sequin_directory" ]; then
	mkdir sequin_directory
fi

# Getting info for virus nomenclature and divergence 
for feat_tbl2 in *.NT.tbl ; do 
	file_core=${feat_tbl2%.NT.tbl}
	echo $file_core
	file_numbers=$( echo ${file_core: -3} | sed 's/[a-z]//g' | sed 's/[A-Z]//g' )
	echo $file_numbers
	tax_info=${feat_tbl2%.NT.tbl}.tax_guide.blastx.out
	echo $tax_info
### There are only a limited number of taxons here. If you want to add additional taxons, please do so after the 'Anelloviridae' line.
	if grep -q "Anellovir" $tax_info ; then
		vir_name=Anelloviridae ;
	elif grep -q "Microvir" $tax_info ; then
		vir_name=Microviridae ;
	elif grep -q "microphage" $tax_info ; then
		vir_name=Microviridae ;
	elif grep -q "uncultured marine virus" $tax_info ; then
		vir_name="Virus" ;
	elif grep -q "Inovir" $tax_info ; then
		vir_name=Inoviridae ;
	elif grep -q "Siphovir" $tax_info ; then
		vir_name=Siphoviridae ;
	elif grep -q "Myovir" $tax_info ; then
		vir_name=Myoviridae ;		
	elif grep -q "unclassified dsDNA phage" $tax_info ; then
		vir_name="Phage" ;
	elif grep -q "unclassified ssDNA virus" $tax_info ; then
		vir_name="CRESS virus" ;
	elif grep -q "Lake Sarah" $tax_info ; then
		vir_name="CRESS virus" ;
	elif grep -q "Avon-Heathcote" $tax_info ; then
		vir_name="CRESS virus" ;
	elif grep -q "Circovir" $tax_info ; then
		vir_name=Circoviridae ;
	elif grep -q "Genomovir" $tax_info ; then
		vir_name=Genomoviridae ;
	elif grep -q "Geminivir" $tax_info ; then
		vir_name=Geminiviridae ;
	elif grep -q "Polyoma" $tax_info ; then
		vir_name=Polyomaviridae ;
	elif grep -q "Papillomavir" $tax_info ; then
		vir_name=Papillomaviridae ;
	elif grep -q "No homologues found" $tax_info ; then
		if  [[ $1 = "-given_linear" ]]; then
			vir_name="genetic element" ;
		else
			vir_name="circular genetic element" ;
		fi
	elif grep -q "Podovir" $tax_info ; then
		vir_name=Podoviridae ;
	elif grep -q "Parvovir" $tax_info ; then
		vir_name=Parvoviridae ;
	elif grep -q "Bacilladnavir" $tax_info ; then
		vir_name=Bacilladnaviridae ;
	elif grep -q "Caudovir" $tax_info ; then
		vir_name=Caudovirales ;
	elif grep -q "phage" $tax_info ; then
		vir_name="Phage" ;
	elif grep -q "plasmid" $tax_info ; then
		vir_name="metagenomic plasmid" ;
	elif grep -q "Bacteria" $tax_info ; then
		vir_name="Phage" ;
	elif grep -q "virus" $tax_info ; then
		vir_name="Virus" ;
	else
		if  [[ $1 = "-given_linear" ]]; then
			vir_name="unclassified element" ;
		else
			vir_name="Circular genetic element" ;
		fi
	fi
	echo $vir_name ;
#	feature_head=$( echo ">Feature" $3"-associated" $vir_name $file_numbers" Table1" )
	fsa_head=$( echo $vir_name " sp." )
	tax_guess=$( tail -n1 ${feat_tbl2%.NT.tbl}.tax_guide.blastx.out ) ; 
	perc_id=$( head -n1 ${feat_tbl2%.NT.tbl}.tax_guide.blastx.out | sed 's/ /-/g' | awk '{FS="\t"; OFS="\t"} {print $2" "$3}' ) ;
	rand_id=$( echo $RANDOM | tr '[0-9]' '[a-zA-Z]' | cut -c 1-2 )

# Editing and transferring tbl file and fasta (fsa) files to sequin directory
if [ -s ${feat_tbl2%.NT.tbl}.HH_BLAST.tbl ]; then 
	echo "$(tput setaf 5)tbl file made from BLAST and HHBLITS output: "${feat_tbl2%.NT.tbl}.HH_BLAST.tbl" will be used for sqn generation$(tput sgr 0)" ; 
	if grep -q "(plasmid)" ${feat_tbl2%.NT.tbl}.tax_guide.blastx.out ; then
		cp ${feat_tbl2%.NT.tbl}.HH_BLAST.tbl sequin_directory/${feat_tbl2%.NT.tbl}.PLASMID.tbl ; 
		bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [organism=" headername "] [mol_type=genomic DNA][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=circular] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "]" ; print $seq }' ${feat_tbl2%.NT.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.NT.tbl}.PLASMID.fsa ;

	else
		cp ${feat_tbl2%.NT.tbl}.HH_BLAST.tbl sequin_directory/${feat_tbl2%.NT.tbl}.tbl ; 
		if  [[ $1 = "-given_linear" ]]; then
			bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [organism=" headername "] [mol_type=genomic DNA][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=linear] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "]" ; print $seq }' ${feat_tbl2%.NT.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.NT.tbl}.fsa ;	
		else
			bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [organism=" headername "] [mol_type=genomic DNA][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=circular] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "]" ; print $seq }' ${feat_tbl2%.NT.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.NT.tbl}.fsa ;	
		fi
	fi

elif [ -s ${feat_tbl2%.NT.tbl}.BLASTP_final.tbl ]; then 
	echo "$(tput setaf 5)tbl file made from BLASTP + RPS-BLAST output only: "${feat_tbl2%.NT.tbl}.BLASTP_final.tbl" will be used for sqn generation.$(tput sgr 0)"
	if grep -q "(plasmid)" ${feat_tbl2%.NT.tbl}.tax_guide.blastx.out ; then
		cp ${feat_tbl2%.NT.tbl}.BLASTP_final.tbl sequin_directory/${feat_tbl2%.NT.tbl}.PLASMID.tbl ; 
		bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [organism=" headername "] [mol_type=genomic DNA][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=circular] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "]" ; print $seq }' ${feat_tbl2%.NT.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.NT.tbl}.PLASMID.fsa ;	
	else
		cp ${feat_tbl2%.NT.tbl}.BLASTP_final.tbl sequin_directory/${feat_tbl2%.NT.tbl}.tbl ; 
		if  [[ $1 = "-given_linear" ]]; then
			bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [organism=" headername "] [mol_type=genomic DNA][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=linear] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "]" ; print $seq }' ${feat_tbl2%.NT.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.NT.tbl}.fsa ;	
		else
			bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [organism=" headername "] [mol_type=genomic DNA][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=circular] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "]" ; print $seq }' ${feat_tbl2%.NT.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.NT.tbl}.fsa ;	
		fi	
	fi
else
	echo "$(tput setaf 5)tbl file made from RPS-BLAST output only: "$feat_tbl2" will be used for sqn generation.$(tput sgr 0)"
	if grep -q "(plasmid)" ${feat_tbl2%.NT.tbl}.tax_guide.blastx.out ; then
		cp ${feat_tbl2%.NT.tbl}.NT_final.tbl sequin_directory/${feat_tbl2%.NT.tbl}.PLASMID.tbl ; 
		bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [organism=" headername "] [mol_type=genomic DNA][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=circular] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "]" ; print $seq }' ${feat_tbl2%.NT.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.NT.tbl}.PLASMID.fsa ;	
	else
		cp ${feat_tbl2%.NT.tbl}.NT_final.tbl sequin_directory/${feat_tbl2%.NT.tbl}.tbl ; 
		if  [[ $1 = "-given_linear" ]]; then
			bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [organism=" headername "] [mol_type=genomic DNA][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=linear] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "]" ; print $seq }' ${feat_tbl2%.NT.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.NT.tbl}.fsa ;	
		else
			bioawk -v srr_var="$srr_number" -v tax_var="$tax_guess" -v perc_var="$perc_id" -v headername="$fsa_head" -v newname="$file_core" -v source_var="$isolation_source" -v rand_var="$rand_id" -v number_var="$file_numbers" -v date_var="$collection_date" -v metgenome_type_var="$metagenome_type" -v srx_var="$srx_number" -v prjn_var="$bioproject" -v samn_var="$biosample" -c fastx '{ print ">" newname " [note= closest relative: " tax_var " " perc_var "] [organism=" headername "] [mol_type=genomic DNA][isolation_source=" source_var "] [isolate=ct" rand_var number_var " ] [country=USA] [collection_date=" date_var "] [metagenome_source=" metgenome_type_var "] [note=genome binned from sequencing reads available in " srx_var "] [topology=circular] [Bioproject=" prjn_var "] [Biosample=" samn_var "] [SRA=" srr_var "]" ; print $seq }' ${feat_tbl2%.NT.tbl}.rotate.fasta > sequin_directory/${feat_tbl2%.NT.tbl}.fsa ;	
		fi	
	fi	
fi ; 
done

#making cmt file for assembly data
for nucl_fa in $novel_fastas ; do
### this feature will only correctly extract coverage info if SPAdes raw outputs are used for input
	coverage=$( head -n1 $nucl_fa | cut -d " " -f 2 | cut -d "_" -f 6 | sed 's/|.*//g' ) 
	echo "StructuredCommentPrefix	##Genome-Assembly-Data-START##" > sequin_directory/${nucl_fa%.fasta}.cmt ;
	echo "Assembly Method	SPAdes v. 3.11.0" >> sequin_directory/${nucl_fa%.fasta}.cmt ;
	echo "Genome Coverage	"$coverage"x" >> sequin_directory/${nucl_fa%.fasta}.cmt ;
	echo "Sequencing Technology	Illumina" >> sequin_directory/${nucl_fa%.fasta}.cmt ;
done


# Running sequin to generate sqn, gbf, and val files for each genome
for nucl_fa in $novel_fastas ; do 
if [ -s "sequin_directory/${nucl_fa%.fasta}.PLASMID.fsa" ] && [ -s "sequin_directory/${nucl_fa%.fasta}.PLASMID.tbl" ]; then
	echo "$(tput setaf 5)necessary files detected for "$nucl_fa" and attempting to use tbl2asn to make sqn file.$(tput sgr 0)" ; 
	mv  sequin_directory/${nucl_fa%.fasta}.cmt sequin_directory/${nucl_fa%.fasta}.PLASMID.cmt
fi
done
### tbl2asn is available from genbank also.
/data/tiszamj/python/linux64.tbl2asn -V vb -t -t $template_file -X C -p sequin_directory/ ;

if  [[ $1 = "-given_circular" ]] ||  [[ $1 = "-given_linear" ]] ||  [[ $1 = "-needs_rotation" ]]; then
	for fsa_file in sequin_directory/*.fsa ; do
		fsa_name2=$( echo ${fsa_file#sequin_directory/} ) ; 
		fsa_name3=$( echo ${fsa_name2%.fsa} | sed 's/.PLASMID//g' )
		seq_name1=$( head -n1 $fsa_name3.fasta | sed 's/>//g; s/|.*//g' | cut -d " " -f1 )
		sed " 1 s/note= closest relative/note= $seq_name1 ; closest relative/" $fsa_file > $fsa_file.temp
		mv $fsa_file.temp $fsa_file
	done
elif [[ $1 = "-default" ]]; then
	for fsa_file in sequin_directory/*.fsa ; do
		fsa_name2=$( echo ${fsa_file#sequin_directory/} ) ; 
		fsa_name3=$( echo ${fsa_name2%.fsa} | sed 's/.PLASMID//g' )
		seq_name1=$( head -n1 $fsa_name3.fasta | sed 's/>//g; s/|.*//g' | cut -d " " -f2 )
		sed " 1 s/note= closest relative/note= $seq_name1 ; closest relative/" $fsa_file > $fsa_file.temp
		mv $fsa_file.temp $fsa_file
	done
fi

/data/tiszamj/python/linux64.tbl2asn -V vb -t -t $template_file -X C -p sequin_directory/ ;

rm *.all_start_stop.txt *.bad_starts.txt *.comb.tbl *.comb2.tbl *.good_start_orfs.txt *.hypo_start_stop.txt *.nucl_orfs.fa *.remove_hypo.txt *.log *.promer.contigs_with_ends.fa *.promer.promer *.out.hhr *.starting_orf.1.fa *.starting_orf.txt *.used_positions.txt 


echo " "
echo "$(tput setaf 3) >>>>>>CENOTE TAKER HAS FINISHED TAKING CENOTES<<<<<< $(tput sgr 0)"




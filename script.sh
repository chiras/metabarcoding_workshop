#!/bin/sh

# variable definition
wd=/Users/ra39huv/tmp_tut/MolBioDiv-MB-2023/Pollinator_16S

vsearch=~/Documents/Arbeit/Projects/Metabarcoding/Bins/vsearch

# sf=/Applications/SeqFilter-2.1.10/bin/SeqFilter-2.1.10
#ps=~/Downloads/python_scripts

sf=~/Documents/Arbeit/Projects/Metabarcoding/Bins/SeqFilter/bin/SeqFilter
ps=~/Documents/Arbeit/Projects/Metabarcoding/Bins/uchime/python_scripts

threads=5

#extracting files
cd $wd

find . -name '*.gz' -print0 | xargs -0 -I {} -P $threads gunzip -d {}

mkdir -p logs

for f in *_R1*.fastq; do

      r=$(sed -e "s/_R1/_R2/" <<< "$f")
      s=$(cut -d_ -f1 <<< "$f")
      p=$(cut -d_ -f2 <<< "$f")
    	total=$(grep "^@M0" $f | wc -l)


      echo " "
      echo "===================================="
      echo "Processing sample $s"
      echo "(F: $f R: $r)"
      echo "===================================="
      $vsearch --fastq_mergepairs  $f --reverse $r --fastq_minovlen 20 --fastq_maxdiffs 10 --fastqout $s.merge.fq --fastq_eeout --relabel R1+2-$s- --threads $threads  2> logs/vsearch.m.$s.log

      var="$(grep "Merged" logs/vsearch.m.$s.log)"
      echo "$s : $var" | tee -a logs/_merging.log

      $vsearch --fastq_filter $s.merge.fq \
            --fastq_maxee 1 \
            --fastq_minlen 200 \
            --fastq_maxlen 550 \
            --fastq_maxns 0 \
            --fastaout $s.mergefiltered.fa \
            --fasta_width 0 --threads $threads 2> logs/vsearch.mf.$s.log

      var="$(grep "sequences kept" logs/vsearch.mf.$s.log)"
      echo "$s : $var" | tee -a logs/_filter.log

done

  echo " "
  echo "===================================="
  echo "ASV generation and mapping"
  echo "===================================="


  cat *mergefiltered.fa > all.merge.fasta

  # cleanup
  mkdir -p raw
  mkdir -p tmp

  mv *.fastq raw/
  mv *.fq tmp/
  mv *.fa tmp/

  #tar

  $vsearch --derep_fulllength all.merge.fasta \
    --minuniquesize 2 \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.merge.derep.uc \
    --output all.merge.derep.fa --threads $threads 2> logs/_derep.log

  $vsearch --cluster_unoise  all.merge.derep.fa --sizein --sizeout --relabel ASV_ --centroids zotus.merge_chim.fa --threads $threads 2> logs/_unoise.log

  grep "Clusters" logs/_unoise.log
  grep "Singletons" logs/_unoise.log

  $vsearch --sortbysize zotus.merge_chim.fa --output zotus.merge_sorted.fa --threads $threads 2>  logs/_sort.log

  $vsearch --uchime3_denovo zotus.merge_sorted.fa --abskew 16 --nonchimeras zotus.merge.fa --threads $threads 2>  logs/_uchime.log

  grep "Found" logs/_uchime.log

  # renaming ASVs

  # [TODO] make decision whether to take merged or R1 here at some point

  ### create community table
  head all.merge.fasta |  sed "s/^>R1+2-\([a-zA-Z0-9-]*\)\-\([0-9]*\)/>R1+2-\1_\2;barcodelabel=\1;/g" > all.merge.bc.fasta

  $vsearch --usearch_global all.merge.bc.fasta --db zotus.merge.fa --strand plus --id 0.97 --uc map.merge.uc --threads $threads 2> logs/_mapping.log

  grep "Matching" logs/_mapping.log

  python2.7 $ps/uc2otutab.py map.merge.uc > asv_table.merge.txt

##### create taxonomy

  echo " "
  echo "===================================="
  echo "Taxonomic classification"
  echo "===================================="

# VSEARCH

# Bacterial references
unset refDBs
refDBs[1]="../_DBs/16S_rdp_16s_v18.fa"
hieDBs="../_DBs/16S_rdp_16s_v18.fa"

# Plant references
#refDBs[1]="../_DBs/its2.Germany_niedersachsen_species.2020-05-11.tax.fa"
#refDBs[2]="../_DBs/its2.viridiplantae.none.2019-07-16.tax.fa"
#hieDBs="../_DBs/its2.viridiplantae.none.2019-07-16.tax.fa"

threshold=97

touch taxonomy.vsearch
echo ",kingdom,phylum,order,family,genus,species" > taxonomy.vsearch

countdb=0
cp  zotus.merge.fa zotus.direct.$countdb.uc.nohit.fasta
prevDB=$countdb

for db in "${refDBs[@]}"
  do :
    countdb=$((countdb+1))
    echo "\n\n#### Direct VSEARCH Classification level: $countdb";
    $vsearch --usearch_global zotus.direct.$prevDB.uc.nohit.fasta --db $db --id 0.$threshold --uc zotus.direct.$countdb.uc --fastapairs zotus.direct.$countdb.fasta --strand both --threads $threads 2>  logs/_direct.$countdb.log

    grep "^N[[:space:]]" zotus.direct.$countdb.uc | cut -f 9 > zotus.direct.$countdb.uc.nohit
    $sf zotus.merge.fa --ids zotus.direct.$countdb.uc.nohit --out zotus.direct.$countdb.uc.nohit.fasta
    cut -f 9,10 zotus.direct.$countdb.uc  | grep -v "*" | sed "s/[A-Za-z0-9_.-]*;tax=//" >> taxonomy.vsearch
    prevDB=$countdb
  done

echo "\n\n#### Hierarchical VSEARCH classification";

$vsearch --sintax zotus.direct.$countdb.uc.nohit.fasta --db $hieDBs --tabbedout zotus.uc.merge.nohit.sintax --strand plus --sintax_cutoff 0.9 --threads $threads 2>  logs/_sintax.log

# for bacteria
cut -f1,4 zotus.uc.merge.nohit.sintax |  sed -E -e "s/,s:.*$//"  >> taxonomy.vsearch

# for plants
# cut -f1,4 zotus.uc.merge.nohit.sintax |  sed -E -e "s/\_[0-9]+//g" -e "s/,s:.*$//"  >> taxonomy.vsearch

sed -i .bak -e "s/c:.*,o:/o:/g" -e "s/[A-Za-z0-9]*;tax=//" -e "s/[[:space:]]/,/" taxonomy.vsearch

sed -i .bak "s/OTUId//" asv_table.merge.txt

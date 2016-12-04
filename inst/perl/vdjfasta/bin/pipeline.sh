hlaFile="tb-patient-hla.txt"
annotationFile="tb-clone-annotations.txt"

hlaFileNew=`echo "$hlaFile" | tr " " "_"`
annotationFileNew=`echo "$annotationFile" | tr " " "_"`


cdr3File=`echo $annotationFile | sed 's/.txt$//' | sed 's/$/.cdr3.txt/'`
cdr3corename=`echo $cdr3File | sed 's/.txt$//'`

cat $annotationFile | grep -v CDR3b | cut -f 1 | sort | uniq > $cdr3File

../bin/gliph-group-discovery.pl --textfile=$cdr3File --lcminp=0.001 --lcminove=10 > $cdr3File-gliph.results.txt

../bin/gliph-group-scoring.pl --convergence_file=$cdr3corename-convergence-groups.txt --clone_annotations=$annotationFile --hla_file=$hlaFile --motif_pval_file=$cdr3corename-kmer_resample_1000_minp0.001_ove10.txt > ${cdr3corename}-gliph-group-scoring.txt



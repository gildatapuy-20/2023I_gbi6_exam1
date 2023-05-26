echo "Extensión de archivos"
echo "TXT:" `find ../exam1/data/singlecell | grep -w "txt" | wc -l`
echo "PNG:" `find ../exam1/data/singlecell | grep -w "png" | wc -l`
echo "PDF:" `find ../exam1/data/singlecell | grep -w "pdf" | wc -l`
echo "GZ:" `find ../exam1/data/singlecell | grep -w "gz" | wc -l`
echo "BAM:" `find ../exam1/data/singlecell | grep -w "bam" | wc -l`
echo "BAI:" `find ../exam1/data/singlecell | grep -w "bai" | wc -l`
echo "TSV:" `find ../exam1/data/singlecell | grep -w "tsv" | wc -l`
echo "FASTQ:" `find ../exam1/data/singlecell | grep -w "fastq" | wc -l`

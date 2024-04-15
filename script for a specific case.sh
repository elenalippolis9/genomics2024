mkdir case600
cd case600




# Align reads and create BAM files for each individual 
bowtie2 -U /home/BCG2024_genomics_exam/case600_father.fq.gz -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SF' --rg "SM:father" | samtools view -Sb | samtools sort -o case600_father.bam
bowtie2 -U /home/BCG2024_genomics_exam/case600_child.fq.gz -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SC' --rg "SM:child" | samtools view -Sb | samtools sort -o case600_child.bam
bowtie2 -U /home/BCG2024_genomics_exam/case600_mother.fq.gz -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SM' --rg "SM:mother" | samtools view -Sb | samtools sort -o case600_mother.bam

# Index BAM files
samtools index case600_father.bam
samtools index case600_child.bam
samtools index case600_mother.bam

#Run fastQC 
fastqc /home/BCG2024_genomics_exam/case600_father.fq.gz
fastqc /home/BCG2024_genomics_exam/case600_child.fq.gz
fastqc /home/BCG2024_genomics_exam/case600_mother.fq.gz

# Run QC using qualimap and rename the stats file
qualimap bamqc -bam case600_child.bam -gff /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -outdir case600_child
#mv "${case_name}_child_stats" "${case_name}_child"
qualimap bamqc -bam case600_mother.bam -gff /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -outdir case600_mother
#mv "${case_name}_mother_stats" "${case_name}_mother"
qualimap bamqc -bam case600_father.bam -gff /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -outdir case600_father
#mv "${case_name}_father_stats" "${case_name}_father"


# Performing multiqc
multiqc .

# Run FreeBayes to generate VCF
nohup freebayes -f /home/BCG2024_genomics_exam/universe.fasta -m 20 -C 5 -Q 10 --min-coverage 10 case600_child.bam case600_father.bam case600_mother.bam > case600.vcf &
wait 

# Order the patients in the VCF file to be sure they match the order of the pattern we put for the variants selection
bcftools query -l case600.vcf | sort > case600.samples.txt
bcftools view -S case600.samples.txt case600.vcf > case600.sorted.vcf

# Show and count the variants of each individual to understand how setting up the segregation patterns
grep -v "#"  case600.sorted.vcf | cut -f 10 | cut -d ":" -f 1 | sort | uniq -c
grep -v "#"  case600.sorted.vcf | cut -f 11 | cut -d ":" -f 1 | sort | uniq -c
grep -v "#"  case600.sorted.vcf | cut -f 12 | cut -d ":" -f 1 | sort | uniq -c

# Select variants based on case type
if [ "$case_type" == "AR" ]; then
    # Select variants with the recessive pattern
    grep "^#" case600.sorted.vcf > candilist600.vcf
    grep "1/1.*0/1.*0/1" case600.sorted.vcf >> candilist600.vcf
elif [ "$case_type" == "AD" ]; then
    # Select variants with the dominant pattern
    grep "^#" case600.sorted.vcf > candilist600.vcf
    grep -E "(0/1.*0/0.*0/0|0/1.*0/1.*0/0|0/1.*0/0.*0/1)" case600.sorted.vcf >> candilist600.vcf
else
    # In case the type (recessive or dominant) is not specified
    echo "Invalid case type: $case_type"
fi

# Generation of the VCF file keeping just the variants included in the target regions (-u serve per riportare le cose una voltaq
# sola in caso di overlapping exons)
grep "^#" candilist600.vcf > 600candilistTG.vcf
bedtools intersect -a candilist600.vcf -b /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -u >> 600candilistTG.vcf

#Generate coverage tracks with bedtools genomecov
bedtools genomecov -ibam case600_father.bam -bg -trackline -trackopts 'name="father"' -max 100 > 600fatherCov.bg
bedtools genomecov -ibam case600_mother.bam -bg -trackline -trackopts 'name="mother"' -max 100 > 600motherCov.bg
bedtools genomecov -ibam case600_child.bam -bg -trackline -trackopts 'name="child"' -max 100 > 600childCov.bg

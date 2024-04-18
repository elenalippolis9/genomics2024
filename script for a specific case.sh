mkdir case644
cd case644

# Align reads and create BAM files for each individual 
bowtie2 -U /home/BCG2024_genomics_exam/case644_father.fq.gz -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SF' --rg "SM:father" | samtools view -Sb | samtools sort -o case644_father.bam
bowtie2 -U /home/BCG2024_genomics_exam/case644_child.fq.gz -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SC' --rg "SM:child" | samtools view -Sb | samtools sort -o case644_child.bam
bowtie2 -U /home/BCG2024_genomics_exam/case644_mother.fq.gz -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SM' --rg "SM:mother" | samtools view -Sb | samtools sort -o case644_mother.bam

# Index BAM files
samtools index case644_father.bam
samtools index case644_child.bam
samtools index case644_mother.bam

# Run fastQC 
mkdir fastqc_outputs
cd fastqc_outputs
fastqc /home/BCG2024_genomics_exam/case644_father.fq.gz -o ./
fastqc /home/BCG2024_genomics_exam/case644_child.fq.gz  -o ./
fastqc /home/BCG2024_genomics_exam/case644_mother.fq.gz -o ./
cd ..

# Run QC using qualimap and rename the stats file
mkdir qualimap_outputs
cd qualimap_outputs
qualimap bamqc -bam /home/BCG_2024_llepore/genomics2024/final_project/case600/case644_child.bam -gff /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -outdir case644_child
qualimap bamqc -bam /home/BCG_2024_llepore/genomics2024/final_project/case600/case644_mother.bam -gff /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -outdir case644_mother
qualimap bamqc -bam /home/BCG_2024_llepore/genomics2024/final_project/case600/case644_father.bam -gff /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -outdir case644_father
cd ..

# Performing multiqc
multiqc . -o multiqc_reports

# Run FreeBayes to generate VCF
nohup freebayes -f /home/BCG2024_genomics_exam/universe.fasta -m 20 -C 5 -Q 10 --min-coverage 10 case644_child.bam case644_father.bam case644_mother.bam > case644.vcf &
wait 

# Order the patients in the VCF file to be sure they match the order of the pattern we put for the variants selection
bcftools query -l case644.vcf | sort > case644.samples.txt
bcftools view -S case644.samples.txt case644.vcf > case644.sorted.vcf

# Show and count the variants of each individual to understand how setting up the segregation patterns
grep -v "#"  case644.sorted.vcf | cut -f 10 | cut -d ":" -f 1 | sort | uniq -c
grep -v "#"  case644.sorted.vcf | cut -f 11 | cut -d ":" -f 1 | sort | uniq -c
grep -v "#"  case644.sorted.vcf | cut -f 12 | cut -d ":" -f 1 | sort | uniq -c

# Select variants based on case type
if [ "$case_type" == "AR" ]; then
    # Select variants with the recessive pattern
    grep "^#" case644.sorted.vcf > candilist644.vcf
    grep "1/1.*0/1.*0/1" case644.sorted.vcf >> candilist644.vcf
elif [ "$case_type" == "AD" ]; then
    # Select variants with the dominant pattern
    grep "^#" case644.sorted.vcf > candilist644.vcf
    grep "0/1.*0/0.*0/0" case644.sorted.vcf >> candilist644.vcf
else
    # In case the type (recessive or dominant) is not specified
    echo "Invalid case type: $case_type"
fi

# Generation of the VCF file keeping just the variants included in the target regions (-u serve per riportare le cose una voltaq
# sola in caso di overlapping exons)
grep "^#" candilist644.vcf > 644candilistTG.vcf
bedtools intersect -a candilist644.vcf -b /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -u >> 644candilistTG.vcf

#Generate coverage tracks with bedtools genomecov
bedtools genomecov -ibam case644_father.bam -bg -trackline -trackopts 'name="father"' -max 100 > 644fatherCov.bg
bedtools genomecov -ibam case644_mother.bam -bg -trackline -trackopts 'name="mother"' -max 100 > 644motherCov.bg
bedtools genomecov -ibam case644_child.bam -bg -trackline -trackopts 'name="child"' -max 100 > 644childCov.bg

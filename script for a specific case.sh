mkdir case717
cd case717

# Align reads and create BAM files for each individual 
bowtie2 -U /home/BCG2024_genomics_exam/case717_father.fq.gz -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SF' --rg "SM:father" | samtools view -Sb | samtools sort -o case717_father.bam
bowtie2 -U /home/BCG2024_genomics_exam/case717_child.fq.gz -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SC' --rg "SM:child" | samtools view -Sb | samtools sort -o case717_child.bam
bowtie2 -U /home/BCG2024_genomics_exam/case717_mother.fq.gz -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SM' --rg "SM:mother" | samtools view -Sb | samtools sort -o case717_mother.bam

# Index BAM files
samtools index case717_father.bam
samtools index case717_child.bam
samtools index case717_mother.bam

# Run fastQC 
mkdir fastqc_outputs
cd fastqc_outputs
fastqc /home/BCG2024_genomics_exam/case717_father.fq.gz -o ./
fastqc /home/BCG2024_genomics_exam/case717_child.fq.gz  -o ./
fastqc /home/BCG2024_genomics_exam/case717_mother.fq.gz -o ./
cd ..

# Run QC using qualimap and rename the stats file
mkdir qualimap_outputs
cd qualimap_outputs
qualimap bamqc -bam /home/BCG_2024_llepore/genomics2024/final_project/case600/case717_child.bam -gff /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -outdir case717_child
qualimap bamqc -bam /home/BCG_2024_llepore/genomics2024/final_project/case600/case717_mother.bam -gff /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -outdir case717_mother
qualimap bamqc -bam /home/BCG_2024_llepore/genomics2024/final_project/case600/case717_father.bam -gff /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -outdir case717_father
cd ..

# Performing multiqc
multiqc . -o multiqc_reports

# Run FreeBayes to generate VCF
nohup freebayes -f /home/BCG2024_genomics_exam/universe.fasta -m 20 -C 5 -Q 10 --min-coverage 10 case717_child.bam case717_father.bam case717_mother.bam > case717.vcf &
wait 

# Order the patients in the VCF file to be sure they match the order of the pattern we put for the variants selection
bcftools query -l case717.vcf | sort > case717.samples.txt
bcftools view -S case717.samples.txt case717.vcf > case717.sorted.vcf

# Show and count the variants of each individual to understand how setting up the segregation patterns
grep -v "#"  case717.sorted.vcf | cut -f 10 | cut -d ":" -f 1 | sort | uniq -c
grep -v "#"  case717.sorted.vcf | cut -f 11 | cut -d ":" -f 1 | sort | uniq -c
grep -v "#"  case717.sorted.vcf | cut -f 12 | cut -d ":" -f 1 | sort | uniq -c

# Select variants based on case type
if [ "$case_type" == "AR" ]; then
    # Select variants with the recessive pattern
    grep "^#" case717.sorted.vcf > candilist717.vcf
    grep "1/1.*0/1.*0/1" case717.sorted.vcf >> candilist717.vcf
elif [ "$case_type" == "AD" ]; then
    # Select variants with the dominant pattern
    grep "^#" case717.sorted.vcf > candilist717.vcf
    grep "0/1.*0/0.*0/0" case717.sorted.vcf >> candilist717.vcf
else
    # In case the type (recessive or dominant) is not specified
    echo "Invalid case type: $case_type"
fi

# Generation of the VCF file keeping just the variants included in the target regions (-u serve per riportare le cose una voltaq
# sola in caso di overlapping exons)
grep "^#" candilist717.vcf > 717candilistTG.vcf
bedtools intersect -a candilist717.vcf -b /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -u >> 717candilistTG.vcf

#Generate coverage tracks with bedtools genomecov
bedtools genomecov -ibam case717_father.bam -bg -trackline -trackopts 'name="father"' -max 100 > 717fatherCov.bg
bedtools genomecov -ibam case717_mother.bam -bg -trackline -trackopts 'name="mother"' -max 100 > 717motherCov.bg
bedtools genomecov -ibam case717_child.bam -bg -trackline -trackopts 'name="child"' -max 100 > 717childCov.bg

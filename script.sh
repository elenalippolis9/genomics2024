cases=("case622 AR" "case717 AD" "case644 AD" "case625 AD" "case743 AD" "case584 AD" "case709 AR" "case600 AR" "case696 AD" "case630 AD")
for case_info in "${cases[@]}"; do
    echo "$case_info"
    #if directory e file già creati move on
    #else create the directory and file with results
    case_name="${case_info%% *}"  # Extracts the substring before the first space
    case_type="${case_info#* }"   # Extracts the substring after the first space
    
   # Align reads and create BAM files for each individual
    bowtie2 -U "$case_name"_father.fq.gz -p 8 -x uni --rg-id 'SF' --rg "SM:father" | samtools view -Sb | samtools sort -o "working/$case_name/case${case_name}_father.bam"
    bowtie2 -U "$case_name"_child.fq -p 8 -x uni --rg-id 'SC' --rg "SM:child" | samtools view -Sb | samtools sort -o "working/$case_name/case${case_name}_child.bam"
    bowtie2 -U "$case_name"_mother.fq.gz -p 8 -x uni --rg-id 'SM' --rg "SM:mother" | samtools view -Sb | samtools sort -o "working/$case_name/case${case_name}_mother.bam"

    # Index BAM files
    samtools index "${case_name}_father.bam"
    samtools index "${case_name}_child.bam"
    samtools index "${case_name}_mother.bam"

    # Run QC using qualimap
    qualimap bamqc -bam "${case_name}_child.bam"
    qualimap bamqc -bam "${case_name}_mother.bam"
    qualimap bamqc -bam "${case_name}_father.bam"

    ./multiqc

    # Generate coverage tracks with bedtools genomecov
    bedtools genomecov -ibam "${case_name}_father.bam" -bg -trackline -trackopts 'name="father"' -max 100 > "father${case_name}Cov.bg"
    bedtools genomecov -ibam "${case_name}_mother.bam" -bg -trackline -trackopts 'name="mother"' -max 100 > "mother${case_name}Cov.bg"
    bedtools genomecov -ibam "${case_name}_child.bam" -bg -trackline -trackopts 'name="child"' -max 100 > "child${case_name}Cov.bg"

    # Run FreeBayes to generate VCF
    nohup freebayes -f universe.fasta -m 20 -C 5 -Q 10 --min-coverage 10 "${case_name}_child.bam" "${case_name}_father.bam" "${case_name}_mother.bam" > "vcf_files/${case_name}.vcf" &

    # Select variants based on case type
    if [ "$case_type" == "AR" ]; then
        # Select variants with the "right" segregation pattern
        grep "^#" "${case_name}.vcf" > "candilist${case_name}.vcf"
        grep "0/1.*0/1.*1/1" "${case_name}.vcf" >> "candilist${case_name}.vcf"

        # Keep only the variants included in the target regions
        bedtools intersect -a "candilist${case_name}.vcf" -b exons16Padded_sorted.bed -u > "${case_name}candilistTG.vcf"
    elif [ "$case_type" == "AD" ]; then
     # Other operations for AD cases
        echo "Performing other operations for AD cases."
    else
        echo "Invalid case type: $case_type"
    fi

done
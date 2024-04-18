cases=("case622 AR" "case717 AD" "case644 AD" "case625 AD" "case743 AD" "case584 AD" "case709 AR" "case600 AR" "case696 AD" "case630 AD")
for case_info in "${cases[@]}"; do
    case_name="${case_info%% *}"  # Extracts the substring before the first space
    case_type="${case_info#* }"   # Extracts the substring after the first space
    individuals=("mother" "father" "child")
    
    directory="/home/BCG_2024_elippolis/${case_name}"
    
    #if directory e file già creati move on
    if [ -d "$directory" ]; then
        echo "Directory and file already exist. Moving on..."

    #else create the directory and file with results
    else
    # Create the directory if it doesn't exist
        mkdir -p "$directory"
    fi

 #nano Align reads and create BAM files for each individual 
  
        bowtie2 -U "/home/BCG2024_genomics_exam/${case_name}_mother.fq.gz" -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SM' --rg "SM:mother" | samtools view -Sb | samtools sort -o "${directory}/${case_name}_mother.bam"
        bowtie2 -U "/home/BCG2024_genomics_exam/${case_name}_father.fq.gz" -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SF' --rg "SM:father" | samtools view -Sb | samtools sort -o "${directory}/${case_name}_father.bam"
        bowtie2 -U "/home/BCG2024_genomics_exam/${case_name}_child.fq.gz" -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SC' --rg "SM:child" | samtools view -Sb | samtools sort -o "${directory}/${case_name}_child.bam"
   
   
    # Index BAM files
    for individual in "${individuals[@]}"; do
   	  samtools index "${directory}/${case_name}_${individual}.bam"
    done
   

    #Run fastQC 
    for individual in "${individuals[@]}"; do
   	  fastqc "/home/BCG2024_genomics_exam/${case_name}_${individual}.fq.gz" -o "${directory}"
    done
    
   
    # Run QC using qualimap and rename the stats file
    for individual in "${individuals[@]}"; do
          qualimap bamqc -bam "${directory}/${case_name}_${individual}.bam" -gff /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -outdir "${directory}/${case_name}_${individual}"

       #  mv "/home/BCG_2024_elippolis/${case_name}_${individual}_stats" "${directory}/${case_name}_${individual}"
    done
    
   
    # Performing multiqc
    multiqc "${directory}" --filename "${directory}/multiqc_report_${case_name}.html"

 # Run FreeBayes to generate VCF
    nohup freebayes -f /home/BCG2024_genomics_exam/universe.fasta -m 20 -C 5 -Q 10 --min-coverage 10 "${directory}/${case_name}_child.bam" "${directory}/${case_name}_father.bam" "${directory}/${case_name}_mother.bam" > "${directory}/${case_name}.vcf" &
    wait
 
    # Order the patients in the VCF file to be sure they match the order of the pattern we put for the variants selection

    bcftools query -l "${directory}/${case_name}.vcf" | sort > "${directory}/${case_name}.samples.txt"
    bcftools view -S "${directory}/${case_name}.samples.txt" "${directory}/${case_name}.vcf" > "${directory}/${case_name}.sorted.vcf"

  
    # Show and count the variants of each individual to understand how setting up the segregation patterns
    grep -v "#"  "${directory}/${case_name}.sorted.vcf" | cut -f 10 | cut -d ":" -f 1 | sort | uniq -c
    grep -v "#"  "${directory}/${case_name}.sorted.vcf" | cut -f 11 | cut -d ":" -f 1 | sort | uniq -c
    grep -v "#"  "${directory}/${case_name}.sorted.vcf" | cut -f 12 | cut -d ":" -f 1 | sort | uniq -c

    grep "^#" "${directory}/${case_name}.sorted.vcf" > "${directory}/candilist${case_name}.vcf"
    # Select variants based on case type
    if [ "$case_type" == "AR" ]; then
        # Select variants with the recessive pattern

        grep "1/1.*0/1.*0/1" "${directory}/${case_name}.sorted.vcf" >> "${directory}/candilist${case_name}.vcf"
    elif [ "$case_type" == "AD" ]; then
        # Select variants with the dominant pattern
        grep "0/1.*0/0.*0/0" "${directory}/${case_name}.sorted.vcf" >> "${directory}/candilist${case_name}.vcf"
    else
        # In case the type (recessive or dominant) is not specified
        echo "Invalid case type: $case_type"
    fi

    # Generation of the VCF file keeping just the variants included in the target regions (-u serve per riportare le cose una voltaq
    # sola in caso di overlapping exons)
    grep "^#" "${directory}/candilist${case_name}.vcf" > "${directory}/${case_name}candilistTG.vcf"
    bedtools intersect -a "${directory}/candilist${case_name}.vcf" -b /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -u >> "${directory}/${case_name}candilistTG.vcf"

    #Generate coverage tracks with bedtools genomecov
     for individual in "${individuals[@]}"; do
     bedtools genomecov -ibam "${directory}/${case_name}_${individual}.bam" -bg -trackline -trackopts 'name="${individual}"' -max 100 > "${directory}/${case_name}${individual}Cov.bg"
   	
    done

    echo "${case_name} finished."

done

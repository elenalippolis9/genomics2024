
# aggiunti fastqc
# da sistemare parte sulla qualità provando i comandi uno alla volta 
# messe graffe al comando di align + spostate virgolette per uniformare stilisticamente
# aggiunto ordinamento campioni
# aggiunto conta varianti per verificare che ci siano solo 0 e 1
# messo in ordine pattern recessivo
# aggiunto pattern dominante 
# aggiunta intersezione con i target genes
# spostata alla fine generazione file cov.bg
# comando multiqc (io sapevo fosse "multiqc .")
# gestire questione directory (come e dove salvare i vari file)
# gestire i percorsi dei file direttamente dal server del prof (dovrei aver fatto, ma controlliamo)
# spostare la parte di qualità alla fine?
# impostare dei for che iterino per padre, madre, bimbo? naaah, non se lo merita



cases=("case622 AR" "case717 AD" "case644 AD" "case625 AD" "case743 AD" "case584 AD" "case709 AR" "case600 AR" "case696 AD" "case630 AD")
for case_info in "${cases[@]}"; do
    echo "$case_info"
    #if directory e file già creati move on
    #else create the directory and file with results
    case_name="${case_info%% *}"  # Extracts the substring before the first space
    case_type="${case_info#* }"   # Extracts the substring after the first space
    
   # Align reads and create BAM files for each individual 
    bowtie2 -U "/home/BCG2024_genomics_exam/${case_name}_father.fq.gz" -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SF' --rg "SM:father" | samtools view -Sb | samtools sort -o "working/$case_name/case${case_name}_father.bam"
    bowtie2 -U "/home/BCG2024_genomics_exam/${case_name}_child.fq.gz" -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SC' --rg "SM:child" | samtools view -Sb | samtools sort -o "working/$case_name/case${case_name}_child.bam"
    bowtie2 -U "/home/BCG2024_genomics_exam/${case_name}_mother.fq.gz" -p 8 -x /home/BCG2024_genomics_exam/uni --rg-id 'SM' --rg "SM:mother" | samtools view -Sb | samtools sort -o "working/$case_name/case${case_name}_mother.bam"

    # Index BAM files
    samtools index "${case_name}_father.bam"
    samtools index "${case_name}_child.bam"
    samtools index "${case_name}_mother.bam"

    #Run fastQC 
    fastqc "/home/BCG2024_genomics_exam/${case_name}_father.fq.gz"
    fastqc "/home/BCG2024_genomics_exam/${case_name}_child.fq.gz"
    fastqc "/home/BCG2024_genomics_exam/${case_name}_mother.fq.gz"

    # Run QC using qualimap and rename the stats file
    qualimap bamqc -bam "${case_name}_child.bam" -gff /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -outdir "${case_name}_child"
    mv "${case_name}_child_stats" "${case_name}_child"
    qualimap bamqc -bam "${case_name}_mother.bam" -gff /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -outdir "${case_name}_mother"
    mv "${case_name}_mother_stats" "${case_name}_mother"
    qualimap bamqc -bam "${case_name}_father.bam" -gff /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -outdir "${case_name}_father"
    mv "${case_name}_father_stats" "${case_name}_father"

    # Performing multiqc
    multiqc .

    # Run FreeBayes to generate VCF
    nohup freebayes -f /home/BCG2024_genomics_exam/universe.fasta -m 20 -C 5 -Q 10 --min-coverage 10 "${case_name}_child.bam" "${case_name}_father.bam" "${case_name}_mother.bam" > "vcf_files/${case_name}.vcf" &
    wait

    # Order the patients in the VCF file to be sure they match the order of the pattern we put for the variants selection
    bcftools query -l "${case_name}.vcf" | sort > "${case_name}.samples.txt"
    bcftools view -S "${case_name}.samples.txt" "${case_name}.vcf" > "${case_name}.sorted.vcf"

    # Show and count the variants of each individual to understand how setting up the segregation patterns
    grep -v "#"  "${case_name}.sorted.vcf" | cut -f 10 | cut -d ":" -f 1 | sort | uniq -c
    grep -v "#"  "${case_name}.sorted.vcf" | cut -f 11 | cut -d ":" -f 1 | sort | uniq -c
    grep -v "#"  "${case_name}.sorted.vcf" | cut -f 12 | cut -d ":" -f 1 | sort | uniq -c

    # Select variants based on case type
    if [ "$case_type" == "AR" ]; then
        # Select variants with the recessive pattern
        grep "^#" "${case_name}.sorted.vcf" > "candilist${case_name}.vcf"
        grep "1/1.*0/1.*0/1" "${case_name}.sorted.vcf" >> "candilist${case_name}.vcf"
    elif [ "$case_type" == "AD" ]; then
        # Select variants with the dominant pattern
        grep "^#" "${case_name}.vcf" > "candilist${case_name}.vcf"
        grep -E "(0/1.*0/0.*0/0)" "${case_name}.vcf" >> "candilist${case_name}.vcf"
    else
        # In case the type (recessive or dominant) is not specified
        echo "Invalid case type: $case_type"
    fi

    # Generation of the VCF file keeping just the variants included in the target regions (-u serve per riportare le cose una voltaq
    # sola in caso di overlapping exons)
    grep "^#" "candilist${case_name}.vcf" > "${case_name}candilistTG.vcf"
    bedtools intersect -a "candilist${case_name}.vcf" -b /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -u >> "${case_name}candilistTG.vcf"

    #Generate coverage tracks with bedtools genomecov
    bedtools genomecov -ibam "${case_name}_father.bam" -bg -trackline -trackopts 'name="father"' -max 100 > "${case_name}fatherCov.bg"
    bedtools genomecov -ibam "${case_name}_mother.bam" -bg -trackline -trackopts 'name="mother"' -max 100 > "${case_name}motherCov.bg"
    bedtools genomecov -ibam "${case_name}_child.bam" -bg -trackline -trackopts 'name="child"' -max 100 > "${case_name}childCov.bg"


done

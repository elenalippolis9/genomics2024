# simil script: simulo di essere all'interno della cartella del caso specifico, andando a creare al suo interno tre sottocartelle:
# fastqc_outputs, qualimap_outputs e multiqc_reports. Alla fine, bisogna solo scaricare il report dalla sottocartella relativa.
# Terminata l'esecuzione, si è nella cartella del caso e si prosegue normalmente (valutare di spostare questa parte tutta alla fine
# per non spezzare il discorso Variant Prioritization)

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

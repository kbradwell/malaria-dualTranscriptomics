#!/bin/bash

f=$1

/Users/kbradwell/Documents/seqtk/seqtk sample -s100 $f 2500000 > subsamp_${f%.gz}

/Users/kbradwell/Documents/seqtk/seqtk sample -s100 ${f%_R1.fastq.gz}_R2.fastq.gz 2500000 > subsamp_${f%_R1.fastq.gz}_R2.fastq

/Users/kbradwell/Documents/hisat2-2.1.0/hisat2 -p 4 --rna-strandness FR --max-intronlen 5000 --fr -k 10 --score-min L,0,-0.4 -x ./multSPP_Genome -1 subsamp_${f%_R1.fastq.gz}_R2.fastq -2 subsamp_${f%.gz} 2> hisat2log_multSPP_${f%_R1.fastq.gz}.log | samtools view -Sbo hisat2_multSPP_Genome_${f%_R1.fastq.gz}.bam - 

rm subsamp_${f%.gz}

rm subsamp_${f%_R1.fastq.gz}_R2.fastq

samtools view -b -F 268 hisat2_multSPP_Genome_${f%_R1.fastq.gz}.bam > bothR1R2mapped_primary_hisat2_multSPP_Genome_${f%_R1.fastq.gz}.bam

samtools view bothR1R2mapped_primary_hisat2_multSPP_Genome_${f%_R1.fastq.gz}.bam|cut -f 3 -d "	"|sort|uniq -c > bothR1R2mapped_primary_hisat2_multSPP_Genome_${f%_R1.fastq.gz}_chromosome_count.txt

samtools view bothR1R2mapped_primary_hisat2_multSPP_Genome_${f%_R1.fastq.gz}.bam|cut -f 3 -d "	"|cut -f 1 -d "_"|sort|uniq -c > bothR1R2mapped_primary_hisat2_multSPP_Genome_${f%_R1.fastq.gz}_species_count.txt

rm hisat2_multSPP_Genome_${f%_R1.fastq.gz}.bam

rm bothR1R2mapped_primary_hisat2_multSPP_Genome_${f%_R1.fastq.gz}.bam

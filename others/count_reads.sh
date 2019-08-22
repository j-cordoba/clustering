#!/bin/bash

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -q smallnodes.q
#$ -pe snode 1
#$ -m beas

#$ -N counting

for i in YSD DRK UNK MTC

do
        grep -c '@SRR' /media/vol2/home/jcordoba/Transcriptomics/"$i"/raw/download/"$i"*.fastq > download_"$i"
        grep -c '@SRR' /media/vol2/home/jcordoba/Transcriptomics/"$i"/raw/prinseq/"$i"*.fastq > prinseq_"$i"
        grep -c '@SRR' /media/vol2/home/jcordoba/Transcriptomics/"$i"/raw/trimmo/"$i"*.fastq > trimmo_"$i"
        grep -c 'SRR' /media/vol2/home/jcordoba/Transcriptomics/"$i"/raw/conta/dna/*.bck > conta_"$i"
        grep -c '>SRR' /media/vol2/home/jcordoba/Transcriptomics/"$i"/raw/reads/{"$i"?_pus_1.fasta,"$i"P?_2.fasta} > reads_"$i"
done

for i in KLY

do
        grep -c '@ERR' /media/vol2/home/jcordoba/Transcriptomics/"$i"/raw/download/"$i"*.fastq > download_"$i"
        grep -c '@ERR' /media/vol2/home/jcordoba/Transcriptomics/"$i"/raw/prinseq/"$i"*.fastq > prinseq_"$i"
        grep -c '@ERR' /media/vol2/home/jcordoba/Transcriptomics/"$i"/raw/trimmo/"$i"*.fastq > trimmo_"$i"
        grep -c 'ERR' /media/vol2/home/jcordoba/Transcriptomics/"$i"/raw/conta/dna/*.bck > conta_"$i"
        grep -c '>ERR' /media/vol2/home/jcordoba/Transcriptomics/"$i"/raw/reads/{"$i"?_pus_1.fasta,"$i"P?_2.fasta} > reads_"$i"
done

for i in EML

do
        ##grep -c '@HISEQ' /media/vol2/home/jcordoba/Transcriptomics/"$i"/raw/download/"$i"*.fastq > download_"$i"
        grep -c '@HISEQ' /media/vol2/home/jcordoba/Transcriptomics/"$i"/raw/prinseq/"$i"*.fastq > prinseq_"$i"
        grep -c '@HISEQ' /media/vol2/home/jcordoba/Transcriptomics/"$i"/raw/trimmo/"$i"*.fastq > trimmo_"$i"
        grep -c 'HISEQ' /media/vol2/home/jcordoba/Transcriptomics/"$i"/raw/conta/dna/*.bck > conta_"$i"
        grep -c '>HISEQ' /media/vol2/home/jcordoba/Transcriptomics/"$i"/raw/reads/{"$i"?_pus_1.fasta,"$i"P?_2.fasta} > reads_"$i"
done

for i in DRK EML KLY MTC UNK YSD

        do
        cat download_"$i" | cut -d ":" -f2 | awk '{s+=$1}END{print s}' > a
        cat prinseq_"$i" | cut -d ":" -f2 | awk '{s+=$1}END{print s}' > b
        cat trimmo_"$i" | cut -d ":" -f2 | awk '{s+=$1}END{print s}' > c
        cat reads_"$i" | cut -d ":" -f2 | awk '{s+=$1}END{print s}' > d
        cat conta_"$i" | cut -d ":" -f2 | awk '{s+=$1}END{print}' > e
        cat {a,b,c,d,e} | awk -v i="$i" 'BEGIN{print i}; {print}' > "$i"
        rm {a,b,c,d,e}
done

printf '%s\n' 'Sample' 'Download' 'Prinseq' 'Trimmo' 'Reads' 'Conta' > rows

paste {rows,DRK,EML,KLY,MTC,UNK,YSD} > counting

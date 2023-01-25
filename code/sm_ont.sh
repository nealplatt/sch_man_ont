conda activate ont_process
cd sm_nanopore

#data/schisto_1_2_2_2022
#data/schisto_1_2_2_2022_renew
#data/schisto_2_9_2022
#data/schisto_frozen_2_12_22
#data/schisto_pairs_04132022
#data/schisto_04122022/

cd results

for run in 


zcat ../data/schisto_frozen_2_12_22/no_sample/20220212_0947_MC-110254_FAR72524_63199f1c/fastq_pass/FAR72524_pass_71f6a84e_*gz \
    >raw_reads/schisto_frozen_2_12_22_raw_pass.fastq &

zcat ../data/schisto_04122022/no_sample/20220412_1058_MC-110254_FAS55645_d05ab9a3/fastq_pass/FAS55645_pass_462f76f0_*gz \
    ../data/schisto_pairs_04132022/no_sample/20220413_1030_MC-110254_FAS57617_fdc4f280/fastq_pass/FAS57617_pass_f77dcaa6_*gz \
    >raw_reads/schisto_pairs_04132022_raw_pass.fastq &

zcat ../data/schisto_1_2_2_2022/no_sample/20220202_1128_MC-110254_FAR74867_c630c15d/fastq_pass/FAR74867_pass_dc8deee6_*gz \
    ../data/schisto_1_2_2_2022_renew/no_sample/20220203_0745_MC-110254_FAR74867_7da50210/fastq_pass/FAR74867_pass_7489bcee_*gz \
    >raw_reads/schisto_1_2_2_2022_raw_pass.fastq &

zcat ../data/schisto_2_9_2022/no_sample/20220209_1005_MC-110254_FAR72529_9c84f0f0/fastq_pass/FAR72529_pass_172de0e0_*gz \
    >raw_reads/schisto_2_9_2022_raw_pass.fastq &

wait

#make histogram and get n50
conda activate n50
n50 \
    --format tsv \
    $(ls *raw_pass.fastq) \
    >n50_raw_cass.tsv

n50 \
    -x \
    $(ls *filtered.fastq)

#.---------------------------------------------------------------------------------------------------------------------------.
#| File                                  | Seqs      | Total bp       | N50    | min | max     | N75    | N90   | auN        |
#+---------------------------------------+-----------+----------------+--------+-----+---------+--------+-------+------------+
#| schisto_pairs_04132022_raw_pass.fastq | 1,490,815 | 13,652,312,137 | 24,594 | 125 | 202,246 | 11,050 | 4,455 | 240,659.34 |
#| schisto_2_9_2022_raw_pass.fastq       | 1,289,677 |  7,662,701,465 | 15,002 | 106 | 143,450 |  5,133 | 2,411 | 214,373.11 |
#| schisto_1_2_2_2022_raw_pass.fastq     | 1,155,730 |  5,567,226,951 | 10,727 |  98 | 133,048 |  3,959 | 1,997 | 206,583.72 |
#| schisto_frozen_2_12_22_raw_pass.fastq | 3,781,415 |  9,440,739,701 |  5,029 | 106 | 136,281 |  2,034 |   968 | 650,254.76 |
#'---------------------------------------+-----------+----------------+--------+-----+---------+--------+-------+------------'


conda deactivate


mkdir ../nanofilt

cd ../nanofilt
for RUN_FASTQ in $(ls ../raw_reads/*.fastq); do
    RUN_ID=$(basename $RUN_FASTQ _raw_pass.fastq)

    NanoFilt \
        --length 500 \
        --quality 7 \
        --headcrop 25 \
        --tailcrop 25 \
        ${RUN_FASTQ} \
        >${RUN_ID}_filtered.fastq &

done

conda activate n50
n50 \
    --format tsv \
    $(ls *filtered.fastq) \
    >n50_filtered.tsv

#.---------------------------------------------------------------------------------------------------------------------------.
#| File                                  | Seqs      | Total bp       | N50    | min | max     | N75    | N90   | auN        |
#+---------------------------------------+-----------+----------------+--------+-----+---------+--------+-------+------------+
#| schisto_pairs_04132022_filtered.fastq | 1,307,890 | 13,518,169,211 | 24,808 | 500 | 202,196 | 11,291 | 4,665 | 239,107.43 |
#| schisto_2_9_2022_filtered.fastq       | 1,199,277 |  7,566,457,104 | 15,278 | 500 | 143,400 |  5,254 | 2,468 | 211,087.68 |
#| schisto_1_2_2_2022_filtered.fastq     | 1,064,057 |  5,477,677,226 | 10,971 | 500 | 132,998 |  4,056 | 2,050 | 202,835.42 |
#| schisto_frozen_2_12_22_filtered.fastq | 3,034,084 |  8,996,160,947 |  5,393 | 500 | 136,231 |  2,233 | 1,125 | 611,852.55 |
#'---------------------------------------+-----------+----------------+--------+-----+---------+--------+-------+------------'

conda deactivate

REF="/master/nplatt/sm_nanopore/data/genome/Smansoni_v7.fa"

mkdir ../minimap2
cd ../minimap2

for QUERY_FQ in $(ls ../nanofilt/*.fastq); do
    RUN_ID=$(basename $QUERY_FQ _filtered.fastq)
    echo $RUN_ID

    minimap2 \
        -a \
        -o ${RUN_ID}.sam \
        -t 12 \
        -x map-ont \
        ${REF} \
        ${QUERY_FQ}

done

#sort bam file
conda activate samtools
for SAM in $(ls *.sam); do
    RUN_ID=$(basename $SAM .sam)
    echo $RUN_ID

    samtools view -Sb -F4 $SAM | samtools sort >${RUN_ID}_sorted.bam
    samtools index ${RUN_ID}_sorted.bam
done
conda deactivate


#sv call
mkdir ../cutesv
cd ../cutesv
for BAM in $(ls ../minimap2/*_sorted.bam); do
    RUN_ID=$(basename $BAM .bam)

    cuteSV \
        --max_cluster_bias_INS 100 \
        --diff_ratio_merging_INS 0.3 \
        --max_cluster_bias_DEL 100 \
        --diff_ratio_merging_DEL 0.3 \
        --threads 12 \
        --min_support 10 \
        --genotype \
        ${BAM} \
        ${REF} \
        ${RUN_ID}.vcf \
        $(pwd)

done

SM_V7_3:21,465,712-21,485,724
SM_V7_3:21,793,415-21,854,235 - polymorphic TE in a gene
SM_V7_3:23,084,116-23,090,873
SM_V7_3:25,973,209-26,003,618

SM_V7_4:28,882,610-28,912,159 - TEs in the UTRs

conda activate repeatmasker
RepeatMasker \
    -pa 12 \
    -lib ../repeatmodeler/Smansoni_v7-families.fa \
    -gff \
    -s \
    -dir $(pwd) \
    -small \
    -x \
    $REF


SM_V7_ZW:1169085
SM_V7_6:4,216,021-4,246,141

Chr: SM_V7_ZW
Position: 1169085
ID: cuteSV.INS.3323
Reference: T*
Alternate: TCATAAGCCGTAGATGTTGTGGCTTGCTCGACAACCGCTGGTGTGATCAAGAAATGCCATCGGTTTATTCTTTGATGAATGATACAGTGAAATATTACGAACTATCAATTACTTACCAAACGGGATTCAGGCTGGTGCAGTTTTCTGATCTATTCCTCAAGTCAGATTCACTTGTACTTGGAGTTGCTGTTGCTAGTTCGACAACTGCTATGTGATCAAACGTTGCCATCTGTCATTCATTGAATGACATGAAAAAATATTACGAACTATCAACAATTACTACTAAACGATTCAGGCGTCAAAAATTTTCTCCGATTCATTCCTCAAGTCAGATTCATCGTGTCTTGTTACGTTACAGTTCAATAACCAAGTTGGCGGATCAAACGTTGCCATTCGAGCTATTCTGATGAATGATATTATGAAACATAAGACCATCAATTACTACCAAACCAGATTCAGTGCTGTGAAGCCCGATCTTATCTCTGCTCAGACCGCGCTGCCATTCGAGCTGGCTCCCGACAATGACACAGCGGAAATATTACGAACGTATCAACACAACACAAATGATTCAGGTGTTTGTGTAGTTTTTCCATCCACTTCAAACGACGCCATGGCCAGATGTTGTGGCTTGTTCGACAACTGCTGGAAAGATCAAAAACGTTGCCATTCGAGCTATTCTGATGA
Qual: 136.5
Type: INDEL
Is Filtered Out: No

Alleles:
Alternate Alleles: TCATAAGCCGTAGATGTTGTGGCTTGCTCGACAACCGCTGGTGTGATCAAGAAATGCCATCGGTTTATTCTTTGATGAATGATACAGTGAAATATTACGAACTATCAATTACTTACCAAACGGGATTCAGGCTGGTGCAGTTTTCTGATCTATTCCTCAAGTCAGATTCACTTGTACTTGGAGTTGCTGTTGCTAGTTCGACAACTGCTATGTGATCAAACGTTGCCATCTGTCATTCATTGAATGACATGAAAAAATATTACGAACTATCAACAATTACTACTAAACGATTCAGGCGTCAAAAATTTTCTCCGATTCATTCCTCAAGTCAGATTCATCGTGTCTTGTTACGTTACAGTTCAATAACCAAGTTGGCGGATCAAACGTTGCCATTCGAGCTATTCTGATGAATGATATTATGAAACATAAGACCATCAATTACTACCAAACCAGATTCAGTGCTGTGAAGCCCGATCTTATCTCTGCTCAGACCGCGCTGCCATTCGAGCTGGCTCCCGACAATGACACAGCGGAAATATTACGAACGTATCAACACAACACAAATGATTCAGGTGTTTGTGTAGTTTTTCCATCCACTTCAAACGACGCCATGGCCAGATGTTGTGGCTTGTTCGACAACTGCTGGAAAGATCAAAAACGTTGCCATTCGAGCTATTCTGATGA
Allele Frequency: 0.7333

Variant Attributes
Allele Frequency: 0.7333
PRECISE: true
RNAMES: NULL
RE: 22
SVTYPE: INS
END: 1169085
SVLEN: 683
CIPOS: [-117, 117]
CILEN: [-115, 115]



whole gene duplication 
SM_V7_1:38,537,253-38,564,035
SM_V7_2:24,364,039-24,422,052
SM_V7_5:9,869,606-9,877,321

gene deletion
SM_V7_3:46,959,111-46,989,520

exon deletion
SM_V7_4:13,269,621-13,283,867

complex
SM_V7_5:6,362,772-6,395,046
##### cuttag pipeline
### Create conda 
conda create -n cuttag python=3
conda activate cuttag
### Create file
mkdir QC PEAKS TRIM BOWTIE DATA
### install
conda install -y sra-tools fastqc trim-galore star hisat2 bowtie2 subread multiqc bwa samtools
conda install trimmomatic
conda install -c bioconda samtools openssl=1.0
conda install -c bioconda macs2 
conda install -c mvdbeek ucsc_tools
conda install -c bioconda bedtools

### QC
path="/home/zhangquan/Projects/omics/cuttag/DATA/fq/"
address="/home/zhangquan/Projects/omics/cuttag"

fastqc -t 2 -o $address/QC $path/*.fq.gz
multiqc $address/QC/*zip -o $address/QC

# build index
# bowtie2
# Filter reliable SAM (q30)
cd BOWTIE
bowtie2_index="/home/zhangquan/Reference/index/bowtie/mm10/mm10"
path="/home/zhangquan/Projects/omics/cuttag/BOWTIE"
file="/home/zhangquan/Projects/omics/cuttag/DATA/fq"
bowtie2 -p 10 -x $bowtie2_index -X 1000 -1 $file/cDC1-1_FKDL210026203-1a_1.clean.fq.gz -2 $file/cDC1-1_FKDL210026203-1a_2.clean.fq.gz -S $path/cDC11.sam
# 16472680 reads; of these:
  # 16472680 (100.00%) were paired; of these:
    # 650670 (3.95%) aligned concordantly 0 times
    # 13042837 (79.18%) aligned concordantly exactly 1 time
    # 2779173 (16.87%) aligned concordantly >1 times
    # ----
    # 650670 pairs aligned concordantly 0 times; of these:
      # 55268 (8.49%) aligned discordantly 1 time
    # ----
    # 595402 pairs aligned 0 times concordantly or discordantly; of these:
      # 1190804 mates make up the pairs; of these:
        # 792372 (66.54%) aligned 0 times
        # 283361 (23.80%) aligned exactly 1 time
        # 115071 (9.66%) aligned >1 times
# 97.59% overall alignment rate
samtools view -b -f 2 -q 30 -o ./cDC11.pairs.bam ./cDC11.sam

bowtie2 -p 10 -x $bowtie2_index -X 1000 -1 $file/cDC1-2_FKDL210026204-1a_1.clean.fq.gz -2 $file/cDC1-2_FKDL210026204-1a_2.clean.fq.gz -S $path/cDC12.sam
# 17129188 reads; of these:
  # 17129188 (100.00%) were paired; of these:
    # 789232 (4.61%) aligned concordantly 0 times
    # 13506576 (78.85%) aligned concordantly exactly 1 time
    # 2833380 (16.54%) aligned concordantly >1 times
    # ----
    # 789232 pairs aligned concordantly 0 times; of these:
      # 52949 (6.71%) aligned discordantly 1 time
    # ----
    # 736283 pairs aligned 0 times concordantly or discordantly; of these:
      # 1472566 mates make up the pairs; of these:
        # 988098 (67.10%) aligned 0 times
        # 357288 (24.26%) aligned exactly 1 time
        # 127180 (8.64%) aligned >1 times
# 97.12% overall alignment rate
samtools view -b -f 2 -q 30 -o cDC12.pairs.bam cDC12.sam

bowtie2 -p 10 -x $bowtie2_index -X 1000 -1 $file/cDC1-c_FKDL210026202-1a_1.clean.fq.gz -2 $file/cDC1-c_FKDL210026202-1a_2.clean.fq.gz -S $path/cDC1c.sam
# 18146020 reads; of these:
  # 18146020 (100.00%) were paired; of these:
    # 817366 (4.50%) aligned concordantly 0 times
    # 14165255 (78.06%) aligned concordantly exactly 1 time
    # 3163399 (17.43%) aligned concordantly >1 times
    # ----
    # 817366 pairs aligned concordantly 0 times; of these:
      # 85705 (10.49%) aligned discordantly 1 time
    # ----
    # 731661 pairs aligned 0 times concordantly or discordantly; of these:
      # 1463322 mates make up the pairs; of these:
        # 890289 (60.84%) aligned 0 times
        # 408921 (27.94%) aligned exactly 1 time
        # 164112 (11.22%) aligned >1 times
# 97.55% overall alignment rate
samtools view -b -f 2 -q 30 -o cDC1c.pairs.bam cDC1c.sam
#### 排序(bam文件)
bowfile="/home/zhangquan/Projects/omics/cuttag/BOWTIE/bam"
samtools sort -o cDC11.pairs.sort.bam $bowfile/cDC11.pairs.bam
samtools sort -o cDC12.pairs.sort.bam $bowfile/cDC12.pairs.bam
samtools sort -o cDC1c.pairs.sort.bam $bowfile/cDC1c.pairs.bam

#### Library complexity
# preseq
preseqfile="/home/zhangquan/Projects/omics/cuttag/BOWTIE/bam"
preseq lc_extrap -e 1e+8 -P -B -D -v -o cDC11.dat $preseqfile/cDC11.pairs.sort.bam 2>cDC11.log 
preseq lc_extrap -e 1e+8 -P -B -D -v -o cDC12.dat $preseqfile/cDC12.pairs.sort.bam 2>cDC12.log 
preseq lc_extrap -e 1e+8 -P -B -D -v -o cDC1c.dat $preseqfile/cDC1c.pairs.sort.bam 2>cDC1c.log 

#7 peak calling
mkdir bed
cd bed
bedtools bamtobed -i $bowfile/cDC11.pairs.bam > cDC11.final.bed
bedtools bamtobed -i $bowfile/cDC12.pairs.bam > cDC12.final.bed
bedtools bamtobed -i $bowfile/cDC1c.pairs.bam > cDC1c.final.bed

# callpeak
beddir="/home/zhangquan/Projects/omics/cuttag/BOWTIE/bed"
outdir="/home/zhangquan/Projects/omics/cuttag/PEAKS/macs2"
macs2/cDC1c
macs2 callpeak -t $beddir/cDC1c.final.bed -g mm -f BEDPE -n cDC1c  --outdir $outdir -q 0.05 -B --SPMR --keep-dup all 2> cD1c.macs2

macs2/cDC11
macs2 callpeak -t $beddir/cDC11.final.bed -g mm -f BEDPE -n cDC1s  --outdir $outdir -q 0.05 -B --SPMR --keep-dup all 2> cD1s.macs2

macs2/cDC2s
macs2 callpeak -t $beddir/cDC12.final.bed -g mm -f BEDPE -n cDC2s  --outdir $outdir -q 0.05 -B --SPMR --keep-dup all 2> cD2s.macs2

# callpeak P0.05
beddir="/home/zhangquan/Projects/omics/cuttag/nore/BOWTIE/bed"
outdir="/home/zhangquan/Projects/omics/cuttag/nore/PEAKS/macs2_p5"
macs2/cDC11
macs2 callpeak -t $beddir/cDC11.final.bed -g mm -f BEDPE -n cDC1s  --outdir $outdir -p 0.05 -B --SPMR --keep-dup all 2> cD1s_5.macs2

macs2/cDC2s
macs2 callpeak -t $beddir/cDC12.final.bed -g mm -f BEDPE -n cDC2s  --outdir $outdir -p 0.05 -B --SPMR --keep-dup all 2> cD2s_5.macs2

### Prepare the data: The First column: Chromosome; The Second column: Start; The third column: End; The Fourth column: Chain
cDCbed="/home/zhangquan/Projects/omics/cuttag/nore/PEAKS/macs2_p5"
mm10="/home/zhangquan/anaconda3/envs/homer/share/homer/data/genomes"

awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' $cDCbed/cDC11/cDC1s_peaks.narrowPeak > scDC11_homer.bed

awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' $cDCbed/cDC12/cDC2s_peaks.narrowPeak > scDC12_homer.bed

### 发现 motif
outdir="/share/Projects/zhangquan/Projects/omics/Homer/New0916nore/motifcDC11"
findMotifsGenome.pl ./scDC11_homer.bed $mm10/mm10 $outdir -size 200 -len 8,10,12

outdir="/share/Projects/zhangquan/Projects/omics/Homer/New0916nore/motifcDC12"
findMotifsGenome.pl ./scDC12_homer.bed $mm10/mm10 $outdir -size 200 -len 8,10,12


#7 TSS
macs="/home/zhangquan/Projects/omics/cuttag/PEAKS/macs2"
bedSort $macs/cDC1s_treat_pileup.bdg cDC1s_treat_pileup.sort.bdg
samtools view -H $preseqfile/cDC11.pairs.sort.bam |perl -ne 'if(/SN:(\S+)\s+LN:(\d+)/){print "$1\t$2\n"}' >cDC11chrom.size
bedClip -truncate cDC1s_treat_pileup.sort.bdg cDC11chrom.size stdout |perl -ane 'print if($F[1]<$F[2])' >cDC11_treat_pileup.bedGraph
bedGraphToBigWig cDC11_treat_pileup.bedGraph cDC11chrom.size cDC11_treat_pileup.bw

bedSort $macs/cDC2s_treat_pileup.bdg cDC2s_treat_pileup.sort.bdg
samtools view -H $preseqfile/cDC12.pairs.sort.bam |perl -ne 'if(/SN:(\S+)\s+LN:(\d+)/){print "$1\t$2\n"}' >cDC12chrom.size
bedClip -truncate cDC2s_treat_pileup.sort.bdg cDC12chrom.size stdout |perl -ane 'print if($F[1]<$F[2])' >cDC12_treat_pileup.bedGraph
bedGraphToBigWig cDC12_treat_pileup.bedGraph cDC12chrom.size cDC12_treat_pileup.bw

bedSort $macs/cDC1c_treat_pileup.bdg cDC1c_treat_pileup.sort.bdg
samtools view -H $preseqfile/cDC1c.pairs.sort.bam |perl -ne 'if(/SN:(\S+)\s+LN:(\d+)/){print "$1\t$2\n"}' >cDC1cchrom.size
bedClip -truncate cDC1c_treat_pileup.sort.bdg cDC1cchrom.size stdout |perl -ane 'print if($F[1]<$F[2])' >cDC1c_treat_pileup.bedGraph
bedGraphToBigWig cDC1c_treat_pileup.bedGraph cDC1cchrom.size cDC1c_treat_pileup.bw

# deeptools
#gtf2gff
gffread mm10.refGene.gtf -o- > mm10.refGene.gff3
zcat refFlat.txt.gz | awk '{print $3"\t"$5"\t"$5"\t"$2"\t"$1"\t"$4}' > mm10_ucsc_refseq.bed 
bed="/home/zhangquan/Reference/gtf/mm10"
gtf="/home/zhangquan/Reference/gtf/mm10"
tssfile="/home/zhangquan/Projects/omics/cuttag/PEAKS/TSS"
computeMatrix reference-point --referencePoint TSS \
-S $tssfile/cDC1c/cDC1c_treat_pileup.bw $tssfile/cDC11/cDC11_treat_pileup.bw $tssfile/cDC12/cDC12_treat_pileup.bw \
-R $gtf/mm10.refGene.gtf \
-a 3000 \
-b 3000 \
-p 8 \
-o 3mtss3kt.matrix.gz \
--missingDataAsZero \
--skipZeros

plotHeatmap \
--heatmapHeight 16 \
--heatmapWidth 4 \
-m 3mtss3kt.matrix.gz \
-o 3maptss3k.heatmap.pdf \
--colorMap Reds \
--legendLocation none 

##### Use deeptools to draw a signal distribution map near TSS
plotProfile -m 3mtss3kt.matrix.gz \
            -out Zeb1_3mtss3k_Profile.pdf
tssdir="/home/zhangquan/Projects/omics/cuttag/PEAKS/TSS"
# awk '$4> 80 {print}' #$tssdir/cDC1c/cDC1c_treat_pileup.bedGraph > #cDC1c.peak.bed  
cd macs2          
cut -f 1-6 cDC1c_peaks.narrowPeak > macs_cDC1c.bed
cut -f 1-6 cDC1s_peaks.narrowPeak > macs_cDC1s.bed         
cut -f 1-6 cDC2s_peaks.narrowPeak > macs_cDC2s.bed
beddir="/home/zhangquan/Projects/omics/cuttag/PEAKS/macs2/"
tssdir="/home/zhangquan/Projects/omics/cuttag/PEAKS/TSS/"
bed="/home/zhangquan/Reference/gtf/mm10"
gtf="/home/zhangquan/Reference/gtf/mm10"
computeMatrix reference-point --referencePoint TSS \
-S $tssdir/cDC1c/cDC1c_treat_pileup.bw $tssdir/cDC12/cDC12_treat_pileup.bw \
-R $gtf/mm10.refGene.gtf \
-a 3000 \
-b 3000 \
-p 8 \
--skipZeros \
-o cDC1c_cDC2s_3ktss.gz
##### Add --perGroup
plotProfile -m 3mtss3kt.matrix.gz \
            -out Zeb1_3mtss3k_Profile_merge.pdf \
            --perGroup

# The repeatability of peak(IDR)
conda install -c bioconda idr
bedtools intersect -a ./cDC1s_peaks.narrowPeak -b ./cDC2s_peaks.narrowPeak | wc -l

bedtools intersect -a ./cDC1s_peaks.narrowPeak -b ./cDC2s_peaks.narrowPeak > cDC12s_intersect_peaks.bed

bedtools intersect -a ./cDC1s_peaks.narrowPeak -b ./cDC2s_peaks.narrowPeak -f 0.5 -F 0.5 | wc -l

idr \
-s ./cDC1s_peaks.narrowPeak ./cDC2s_peaks.narrowPeak \
-o cDC1s_2s_idr \
--plot


idr \
--samples ./macs2/cDC1s_peaks.narrowPeak ./macs2/cDC2s_peaks.narrowPeak \
--output-file ./cDC1s_cDC2s_idr_narrowPeak.txt \
--plot \
--idr-threshold 0.05

##### Consence Peaks bed file，motif MEME
macsbed="/home/zhangquan/Projects/omics/cuttag/PEAKS/macs2"
mm10fa="/home/zhangquan/Reference/fa"
cut -f 1,2,3 cDC1s_cDC2s_idr_narrowPeak.txt > consence_idr_peaks.bed
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed ./consence_idr_peaks.bed -fo consence_idr_peaks.fa

# There are two ways to handle the results of the output file:
# Through the –soft-idr-threshold parameter to output and save the results of all peaks, the peak of IDR < 0.05 is counted by default, in this case the –output-file parameter saves the .txt file, and the 12th column of the output file records the value of the global IDR value corresponding to peak, and filter by this value;
# The number of peaks output is adjusted by directly defining the IDR value through the –idr-threshold parameter, and the output file contains only the peaks that meet the criteria.

#### Motif MEME AME
# AME is a MEME suite
# Create an environment and install MEME
conda create -n meme
conda install -c conda-forge ucx
conda install -c bioconda perl-xml-parser
conda install -c bioconda meme
conda install -c bioconda bedtools
conda install "python=2.7"
conda install 'icu=58.2'
# motif
# 1. get fa
mm10fa="/home/zhangquan/Reference/fa"
cD1cbed="/home/zhangquan/Projects/omics/cuttag/PEAKS/macs2"
bedtools getfasta -name -fi $mm10fa/mm10.fa -bed $cD1cbed/cDC1c_peaks.narrowPeak -fo cDC1c_peaks.fa
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed $cD1cbed/cDC1s_peaks.narrowPeak -fo cDC11_peaks.fa
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed $cD1cbed/cDC2s_peaks.narrowPeak -fo cDC12_peaks.fa

# Use motif database files for motif analysis
dir="/share/Projects/zhangquan/Projects/omics/"
meme-chip -meme-p 8 -oc $dir/JASresults/ -db ./JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme /home/zhangquan/Projects/omics/cuttag/FA/cDC11_peaks.fa -meme-nmotifs 10

#### Use Homer to find motifs
conda create -n homer
conda install -c bioconda homer
conda install -c bioconda bedtools

### mm10
perl /home/zhangquan/anaconda3/envs/homer/share/homer/configureHomer.pl -install mm10
# downloaded
ls /home/zhangquan/anaconda3/envs/homer/share/homer/data/

### Prepare the data: The First column: Chromosome; The Second column: Start; The third column: End; The Fourth column: Chain
cDCbed="/home/zhangquan/Projects/omics/cuttag/PEAKS/macs2/"
mm10="/home/zhangquan/anaconda3/envs/homer/share/homer/data/genomes/"

awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' $cDCbed/cDC1s_peaks.narrowPeak > cDC11_homer.bed

awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' $cDCbed/cDC2s_peaks.narrowPeak > cDC12_homer.bed

### Find motif
outdir="/share/Projects/zhangquan/Projects/omics/Homer/motif/"
findMotifsGenome.pl ./cDC11_homer.bed $mm10/mm10 $outdir -size 200 -len 8,10,12

outdir="/share/Projects/zhangquan/Projects/omics/Homer/motif12/"
findMotifsGenome.pl ./cDC12_homer.bed $mm10/mm10 $outdir -size 200 -len 8,10,12

##### pvaule & peaks
macsbed="/home/zhangquan/Projects/omics/cuttag/PEAKS/macs2"
mm10fa="/home/zhangquan/Reference/fa"
awk '$7>=3.48768 && $9>=4.17241 {print}' cDC1s_peaks.narrowPeak > cDC1s_subpeaks.narrowPeak.bed

awk '$7>=3.48768 && $9>=4.17241 {print}' cDC2s_peaks.narrowPeak > cDC2s_subpeaks.narrowPeak.bed

bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed $macsbed/cDC1s_subpeaks.narrowPeak.bed -fo cDC11_subpeaks.fa
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed $macsbed/cDC2s_subpeaks.narrowPeak.bed -fo cDC12_subpeaks.fa

##### duplication_peaks_FCmax.txt
macsbed="/home/zhangquan/Projects/omics/cuttag/PEAKS/macs2"
mm10fa="/home/zhangquan/Reference/fa"
sed -i '1d' duplication_peaks_FCmax.txt 
cut -f 1,2,3,6 duplication_peaks_FCmax.txt > dup_peaks_FCmax.bed
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed ./dup_peaks_FCmax.bed -fo dup_peaks_FCmax_intersect.fa

sed -i '1d' duplication_peaks_FCmax_zeb1.txt
cut -f 1,2,3,6 duplication_peaks_FCmax_zeb1.txt > dup_peaks_FCmax_zeb1.bed
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed ./dup_peaks_FCmax_zeb1.bed -fo dup_peaks_FCmax_zeb1_intersect.fa

sed -i '1d' Zeb1c1_peaks.txt
cut -f 1,2,3,4 Zeb1c1_peaks.txt > Zeb1c1.bed
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed ./Zeb1c1.bed -fo Zeb1_cdc1_fj.fa

sed -i '1d' Zeb1c2_peaks.txt
cut -f 1,2,3,4 Zeb1c2_peaks.txt > Zeb1c2.bed
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed ./Zeb1c2.bed -fo Zeb1_cdc2_fj.fa


##### subpeaks
macsbed="/home/zhangquan/Projects/omics/cuttag/PEAKS/macs2"
mm10fa="/home/zhangquan/Reference/fa"
awk '$7>=3.48768 && $9>=4.17241 {print}' cDC1s_peaks.narrowPeak > cDC1s_subpeaks.narrowPeak.txt

awk '$7>=3.48768 && $9>=4.17241 {print}' cDC2s_peaks.narrowPeak > cDC2s_subpeaks.narrowPeak.txt

sed -i '1d' subpeakscDC1_subcDC2_intersect.txt
cut -f 1,2,3,4 subpeakscDC1_subcDC2_intersect.txt > subcDC12_intersect.bed
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed ./subcDC12_intersect.bed -fo subcDC12_inter_fj.fa


##### The overlapping region of two parallel samples
macsbed="/home/zhangquan/Projects/omics/cuttag/PEAKS/macs2"
mm10fa="/home/zhangquan/Reference/fa"
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed $macsbed/cDC12s_intersect_peaks.bed -fo intercDC12_peaks.fa

##### qvaule&FC
macsbed="/home/zhangquan/Projects/omics/cuttag/PEAKS/macs2"
mm10fa="/home/zhangquan/Reference/fa"
sed -i '1d' cdc1_inter_peaks.txt
cut -f 1,2,3,6 cdc1_inter_peaks.txt > cdc1_inter_peaks.bed
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed ./cdc1_inter_peaks.bed -fo cdc1_inter_peaks.fa

sed -i '1d' cdc2_inter_peaks.txt
cut -f 1,2,3,6 cdc2_inter_peaks.txt > cdc2_inter_peaks.bed
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed ./cdc2_inter_peaks.bed -fo cdc2_inter_peaks.fa


sed -i '1d' zeb1_cdc1_pfc_cutoff0427.txt
cut -f 1,2,3,6 zeb1_cdc1_pfc_cutoff0427.txt > cdc1_zeb1_peaks_0427.bed
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed ./cdc1_zeb1_peaks_0427.bed -fo cdc1_zeb1_peaks_0427.fa


sed -i '1d' zeb1_cdc2_pfc_cutoff0427.txt
cut -f 1,2,3,6 zeb1_cdc2_pfc_cutoff0427.txt > cdc2_zeb1_peaks_0427.bed
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed ./cdc2_zeb1_peaks_0427.bed -fo cdc2_zeb1_peaks_0427.fa

sed -i '1d' zeb1_cdc1_fullintersect_0428.txt
cut -f 1,2,3,6 zeb1_cdc1_fullintersect_0428.txt > cdc1_zeb1_peaks_0428.bed
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed ./cdc1_zeb1_peaks_0428.bed -fo cdc1_zeb1_peaks_0428.fa


sed -i '1d' zeb1_cdc2_fullintersect_0428.txt
cut -f 1,2,3,6 zeb1_cdc2_fullintersect_0428.txt > cdc2_zeb1_peaks_0428.bed
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed ./cdc2_zeb1_peaks_0428.bed -fo cdc2_zeb1_peaks_0428.fa

sed -i '1d' cDC1s.peakAnno.txt
cut -f 1,2,3,6 cDC1s.peakAnno.txt > cdc1_full_peaks_0505.bed
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed ./cdc1_full_peaks_0505.bed -fo cdc1_fullanno_peaks_0505.fa

sed -i '1d' zeb1_cdc1_pfc_cutoff0513.txt
cut -f 1,2,3,6 zeb1_cdc1_pfc_cutoff0513.txt > zeb1_cdc1_pfc_cutoff0513.bed
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed ./zeb1_cdc1_pfc_cutoff0513.bed -fo zeb1_cdc1_pfc_cutoff0513.fa

sed -i '1d' zeb1_cdc2_pfc_cutoff0513.txt
cut -f 1,2,3,6 zeb1_cdc2_pfc_cutoff0513.txt > zeb1_cdc2_pfc_cutoff0513.bed
bedtools getfasta  -name -fi $mm10fa/mm10.fa -bed ./zeb1_cdc2_pfc_cutoff0513.bed -fo zeb1_cdc2_pfc_cutoff0513.fa

##### snailpeaks
outdir="/share/Projects/zhangquan/Projects/omics/Homer/New0427snailpeaks/cDC11motif/"
mm10="/home/zhangquan/anaconda3/envs/homer/share/homer/data/genomes"
sed -i '1d' cDC11_snail_peak.txt
cat cDC11_snail_peak.txt > cDC11_snail_homer.bed
findMotifsGenome.pl ./cDC11_snail_homer.bed $mm10/mm10 $outdir -size 200 -len 8,10,12

outdir="/share/Projects/zhangquan/Projects/omics/Homer/New0427snailpeaks/cDC12motif/"
mm10="/home/zhangquan/anaconda3/envs/homer/share/homer/data/genomes"
sed -i '1d' cDC12_snail_peak.txt
cat cDC12_snail_peak.txt > cDC12_snail_homer.bed
findMotifsGenome.pl ./cDC12_snail_homer.bed $mm10/mm10 $outdir -size 200 -len 8,10,12

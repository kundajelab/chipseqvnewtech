#!/usr/bin/env bash

#--------------------
#IMPORTANT: When using .bigWig files, make sure you use the files
# that have "5p" and "counts" in the filename! These are the files corresponding
# to the five-prime ends of the reads. There are also files for coverage, and
# files that are normalized by RPM, but for training models we want the
# *five-prime ends of the raw counts*. Also pay attention to whether the
# .bigWig file has strand-specific counts, and download accordingly (we want
# strand-specific counts for ChIP-seq, but not for CUT-N-RUN)
#--------------------

#We are comparing CUT-N-RUN to ChIP-seq for NANOG and POU5f1 in H1ESCs
#The CUT-N-RUN data is from the paper "Pioneer Factor-Nucleosome Binding Events during Differentiation Are Motif Encoded"

download_utilities () {
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigMerge -O bigWigMerge
    chmod a+x bigWigMerge
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig -O bedGraphToBigWig
    chmod a+x bedGraphToBigWig
}

mergebigwig () {
  local merge1=$1
  local merge2=$2
  local mergeout=$3
  ./bigWigMerge $merge1.bigWig $merge2.bigWig $mergeout.bedGraph
  cat $mergeout.bedGraph | sortBed > $mergeout.sorted.bedGraph
  ./bedGraphToBigWig $mergeout.sorted.bedGraph hg38.chrom.sizes $mergeout.bigWig
  rm $mergeout.bedGraph
  rm $mergeout.sorted.bedGraph
}

bamtobigwig () {
  local inbam=$1
  #get the locations of the 5' reads
  bedtools genomecov -5 -bg -g hg38.chrom.sizes -ibam $inbam.bam | sortBed > $inbam.sorted.bedGraph
  ./bedGraphToBigWig $inbam.sorted.bedGraph hg38.chrom.sizes $inbam.bigWig
  rm $inbam.sorted.bedGraph
}

gather_relevant_files () {
    #Several of the files we will be starting from were prepared by
    # Georgi Marinov; we will just symlink to these.

    local basedir=$1

    #CUT-N-RUN files

    #For the CUT-N-RUN files, e.g. "H1_Nanog-GSM3677838",
    # GSM3677838 is the GEO accession number. maxFL120bp means
    # that the file is filtered for reads where the maximum fragment len is 120bp,
    # and minFL150bp is filtered for reads where the minimum fragment len is 150bp.
    #For CUT-N-RUN, we do NOT want the .bigWig files that separate out the
    # counts by the plus and minus strand
    #(Note: for Pou5f1 in the high Calcium buffer condition are also available,
    #  but they had very few peak calls, so we won't use those)
    #For Pou5f1
    #CUT-N-RUN profiles
    ln -s $basedir/H1_POU5F1_std-GSM3677840/H1_POU5F1_std-GSM3677840.2x25mers.hg38-no_haps.unique.5p.counts.bigWig CUTNRUN.POU5F1.5p.counts.bigWig 
    ln -s $basedir/H1_POU5F1_std-GSM3677840/H1_POU5F1_std-GSM3677840.2x25mers.hg38-no_haps.unique.maxFL120bp.5p.counts.bigWig CUTNRUN.POU5F1.maxFL120bp.5p.counts.bigWig 
    ln -s $basedir/H1_POU5F1_std-GSM3677840/H1_POU5F1_std-GSM3677840.2x25mers.hg38-no_haps.unique.minFL150bp.5p.counts.bigWig CUTNRUN.POU5F1.minFL150bp.5p.counts.bigWig 
    #For the peak file - gzip this one
    cat $basedir/H1_POU5F1_std-GSM3677840/H1_POU5F1_std-GSM3677840.2x25mers.hg38-no_haps.unique.MACS2.IDR-0.05 | gzip -c > CUTNRUN.POU5F1.narrowpeak.gz
    #Get the standard control bigwigs
    #If you look up the GEO accession associated with the POU5f1 experiment, it
    # says that the protocol was "standard" CUT-N-RUN - therefore, we will download
    # the control associated with it.
    ln -s $basedir/H1_IgG_std-GSM3677837/H1_IgG_std-GSM3677837.2x25mers.hg38-no_haps.unique.5p.counts.bigWig CUTNRUN.IgG.standard.5p.counts.bigWig
    ln -s $basedir/H1_IgG_std-GSM3677837/H1_IgG_std-GSM3677837.2x25mers.hg38-no_haps.unique.maxFL120bp.5p.counts.bigWig CUTNRUN.IgG.standard.maxFL120bp.5p.counts.bigWig
    ln -s $basedir/H1_IgG_std-GSM3677837/H1_IgG_std-GSM3677837.2x25mers.hg38-no_haps.unique.minFL150bp.5p.counts.bigWig CUTNRUN.IgG.standard.minFL150bp.5p.counts.bigWig    
    #For NANOG
    #CUT-N-RUN profiles
    ln -s $basedir/H1_Nanog-GSM3677838/H1_Nanog-GSM3677838.2x25mers.hg38-no_haps.unique.5p.counts.bigWig CUTNRUN.NANOG.5p.counts.bigWig
    ln -s $basedir/H1_Nanog-GSM3677838/H1_Nanog-GSM3677838.2x25mers.hg38-no_haps.unique.maxFL120bp.5p.counts.bigWig CUTNRUN.NANOG.maxFL120bp.5p.counts.bigWig
    ln -s $basedir/H1_Nanog-GSM3677838/H1_Nanog-GSM3677838.2x25mers.hg38-no_haps.unique.minFL150bp.5p.counts.bigWig CUTNRUN.NANOG.minFL150bp.5p.counts.bigWig
    #Compress peak file
    cat $basedir/H1_Nanog-GSM3677838/H1_Nanog-GSM3677838.2x25mers.hg38-no_haps.unique.MACS2.IDR-0.05 | gzip -c > CUTNRUN.NANOG.narrowpeak.gz
    #Get both 'auto' control bigwigs
    #If you look up the GEO accession associated with the Nanog experiment, it
    # says that the protocol was "auto" CUT-N-RUN - therefore, we will download
    # both the 'auto' control files and merge them together.
    #IgG A
    ln -s $basedir/H1_IgG_A_auto-GSM3677835/H1_IgG_A_auto-GSM3677835.2x25mers.hg38-no_haps.unique.5p.counts.bigWig CUTNRUN.IgG.A.auto.5p.counts.bigWig
    ln -s $basedir/H1_IgG_A_auto-GSM3677835/H1_IgG_A_auto-GSM3677835.2x25mers.hg38-no_haps.unique.maxFL120bp.5p.counts.bigWig CUTNRUN.IgG.A.auto.maxFL120bp.5p.counts.bigWig
    ln -s $basedir/H1_IgG_A_auto-GSM3677835/H1_IgG_A_auto-GSM3677835.2x25mers.hg38-no_haps.unique.minFL150bp.5p.counts.bigWig CUTNRUN.IgG.A.auto.minFL150bp.5p.counts.bigWig
    #IgG B
    ln -s $basedir/H1_IgG_B_auto-GSM3677836/H1_IgG_B_auto-GSM3677836.2x25mers.hg38-no_haps.unique.5p.counts.bigWig CUTNRUN.IgG.B.auto.5p.counts.bigWig
    ln -s $basedir/H1_IgG_B_auto-GSM3677836/H1_IgG_B_auto-GSM3677836.2x25mers.hg38-no_haps.unique.maxFL120bp.5p.counts.bigWig CUTNRUN.IgG.B.auto.maxFL120bp.5p.counts.bigWig
    ln -s $basedir/H1_IgG_B_auto-GSM3677836/H1_IgG_B_auto-GSM3677836.2x25mers.hg38-no_haps.unique.minFL150bp.5p.counts.bigWig CUTNRUN.IgG.B.auto.minFL150bp.5p.counts.bigWig

    #ChIP-seq files 

    #For ChIP-seq, the files were prepared from the corresponding .bam alignment
    # files, downloaded from the ENCODE portal. For example, in the file name
    # "H1-POU5F1-ENCSR000BMU-ENCFF488NSO.5p.counts.minus.bigWig",
    # ENCSR000BMU is the experiment identifier and ENCFF488NSO is the ID of the .bam
    # file from which the count .bigWig files were prepared.
    #Since the files were prepared by a human under time pressure, you should sanity-check
    # as many things as you can. For example, you for the ChIP-seq files, you can
    # verify that the .bam file corresponding to the file id was aligned to GRCh38.
    # You can also verify that the *filtered* .bam file was used
    # (i.e. "Output Type" is 'alignments' rather than 'unfiltered alignments')
    #For the ChIP-seq, we are interested in the plus and minus strands separately
    #for POU5F1
    ln -s $basedir/ENCODE-H1-ChIP-seq-datasets/H1-POU5F1-ENCSR000BMU-ENCFF488NSO.5p.counts.minus.bigWig ChIPseq.POU5F1.ENCFF488NSO.5p.counts.minus.bigWig
    ln -s $basedir/ENCODE-H1-ChIP-seq-datasets/H1-POU5F1-ENCSR000BMU-ENCFF488NSO.5p.counts.plus.bigWig ChIPseq.POU5F1.ENCFF488NSO.5p.counts.plus.bigWig
    ln -s $basedir/ENCODE-H1-ChIP-seq-datasets/H1-POU5F1-ENCSR000BMU-ENCFF777DLI.5p.counts.minus.bigWig ChIPseq.POU5F1.ENCFF777DLI.5p.counts.minus.bigWig
    ln -s $basedir/ENCODE-H1-ChIP-seq-datasets/H1-POU5F1-ENCSR000BMU-ENCFF777DLI.5p.counts.plus.bigWig ChIPseq.POU5F1.ENCFF777DLI.5p.counts.plus.bigWig
    #We will also get the optimal IDR thresholded peaks
    #The POU5F1 ENCSR000BMU experiment was not put
    # through the ENCODE pipeline to produce an IDR optimal file, so we will use
    # the file Georgi helpfully prepared
    ln -s $basedir/ENCODE-H1-ChIP-seq-datasets/H1-POU5F1-ENCSR000BMU.optimal_idr_thresholded.bed.gz ChIPseq.POU5F1.narrowpeak.gz
    #for NANOG
    ln -s $basedir/ENCODE-H1-ChIP-seq-datasets/H1-NANOG-ENCSR000BMT-ENCFF309DIF.5p.counts.minus.bigWig ChIPseq.NANOG.ENCFF309DIF.5p.counts.minus.bigWig
    ln -s $basedir/ENCODE-H1-ChIP-seq-datasets/H1-NANOG-ENCSR000BMT-ENCFF309DIF.5p.counts.plus.bigWig ChIPseq.NANOG.ENCFF309DIF.5p.counts.plus.bigWig
    ln -s $basedir/ENCODE-H1-ChIP-seq-datasets/H1-NANOG-ENCSR000BMT-ENCFF805SCJ.5p.counts.minus.bigWig ChIPseq.NANOG.ENCFF805SCJ.5p.counts.minus.bigWig
    ln -s $basedir/ENCODE-H1-ChIP-seq-datasets/H1-NANOG-ENCSR000BMT-ENCFF805SCJ.5p.counts.plus.bigWig ChIPseq.NANOG.ENCFF805SCJ.5p.counts.plus.bigWig
    #IDR optimal file for the NANOG ENCSR000BMT experiment
    #For NANOG, the experiment was put through the ENCODE pipeline,
    # so we will just download the IDR optimal file from there
    wget https://www.encodeproject.org/files/ENCFF794GVQ/@@download/ENCFF794GVQ.bed.gz -O ChIPseq.NANOG.narrowpeak.gz
    #For the ChIP-seq control, we will download the .bam file associated with the control
    # for H1ESCs listed on the ENCODE experiment pages for the NANOG and Pou5f1
    # experiments. We will convert the .bam file to a bigwig file.
    wget https://www.encodeproject.org/files/ENCFF948LKH/@@download/ENCFF948LKH.bam -O ChIPseq.control.bam  

    #download the chromsizes file from the DCC
    wget https://raw.githubusercontent.com/ENCODE-DCC/encValData/master/GRCh38/GRCh38_EBV.chrom.sizes -O hg38.chrom.sizes
}

postprocess_files () {
    #We have to mere:
    # the CUT-N-RUN IgG 'auto' controls.
    mergeout=CUTNRUN.IgG.merged.auto.5p.counts
    [[ -f $mergeout.bigWig ]] || mergebigwig CUTNRUN.IgG.A.auto.5p.counts CUTNRUN.IgG.B.auto.5p.counts $mergeout
    mergeout=CUTNRUN.IgG.merged.auto.maxFL120bp.5p.counts
    [[ -f $mergeout.bigWig ]] || mergebigwig CUTNRUN.IgG.A.auto.maxFL120bp.5p.counts CUTNRUN.IgG.B.auto.maxFL120bp.5p.counts $mergeout
    mergeout=CUTNRUN.IgG.merged.auto.minFL150bp.5p.counts
    [[ -f $mergeout.bigWig ]] || mergebigwig CUTNRUN.IgG.A.auto.minFL150bp.5p.counts CUTNRUN.IgG.B.auto.minFL150bp.5p.counts $mergeout
    # the ChIP-seq Pou5f1 replicates
    mergeout=ChIPseq.POU5F1.merged.5p.counts.plus
    [[ -f $mergeout.bigWig ]] || mergebigwig ChIPseq.POU5F1.ENCFF488NSO.5p.counts.plus ChIPseq.POU5F1.ENCFF777DLI.5p.counts.plus $mergeout
    mergeout=ChIPseq.POU5F1.merged.5p.counts.minus
    [[ -f $mergeout.bigWig ]] || mergebigwig ChIPseq.POU5F1.ENCFF488NSO.5p.counts.minus ChIPseq.POU5F1.ENCFF777DLI.5p.counts.minus $mergeout
    # the ChIP-seq NANOG replicates
    mergeout=ChIPseq.NANOG.merged.5p.counts.plus
    [[ -f $mergeout.bigWig ]] || mergebigwig ChIPseq.NANOG.ENCFF309DIF.5p.counts.plus ChIPseq.NANOG.ENCFF805SCJ.5p.counts.plus $mergeout
    mergeout=ChIPseq.NANOG.merged.5p.counts.minus
    [[ -f $mergeout.bigWig ]] || mergebigwig ChIPseq.NANOG.ENCFF309DIF.5p.counts.minus ChIPseq.NANOG.ENCFF805SCJ.5p.counts.minus $mergeout

    #We have to convert the H1 ChIP-seq control from bam to bigwig
    [[ -f ChIPseq.control.bigWig ]] || bamtobigwig ChIPseq.control
}

prepare_bed_regions () {
    #first concatenate all peaks in the CUTNRUN and ChIP-seq files
    zcat CUTNRUN.POU5F1.narrowpeak.gz CUTNRUN.NANOG.narrowpeak.gz ChIPseq.POU5F1.narrowpeak.gz ChIPseq.NANOG.narrowpeak.gz | gzip -c > concatenated_peaks.narrowpeak.gz
    #get 2kb around summits
    zcat concatenated_peaks.narrowpeak.gz | perl -lane 'print $F[0]."\t".(($F[1]+$F[9]))."\t".(($F[1]+$F[9]))' | bedtools slop -g hg38.chrom.sizes -b 1000 | perl -lane 'if ($F[2]-$F[1]==2000) {print $F[0]."\t".$F[1]."\t".$F[2]."\t1"}' | sortBed | gzip -c > 2k_around_summits.bed.gz
    #split into train, valid, test sets by chromosome
    zcat 2k_around_summits.bed.gz | egrep -w 'chr1|chr8|chr21' | gzip -c > test_2k_around_summits.bed.gz
    zcat 2k_around_summits.bed.gz | egrep -w 'chr22' | gzip -c > valid_2k_around_summits.bed.gz
    zcat 2k_around_summits.bed.gz | egrep -w -v 'chr1|chr8|chr21|chr22' | gzip -c > train_2k_around_summits.bed.gz
}

#gather_relevant_files /data/chipseq_vs_newtech/2019-06-26-CUTnRUN-EChO
#download_utilities
#postprocess_files
prepare_bed_regions


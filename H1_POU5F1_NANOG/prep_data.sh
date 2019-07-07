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

#We are comparing CUT-N-RUN to ChIP-seq for Nanog and POU5f1 in H1ESCs
#The CUT-N-RUN data is from the paper "Pioneer Factor-Nucleosome Binding Events during Differentiation Are Motif Encoded"

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
    #For Nanog
    #CUT-N-RUN profiles
    ln -s $basedir/H1_Nanog-GSM3677838/H1_Nanog-GSM3677838.2x25mers.hg38-no_haps.unique.5p.counts.bigWig CUTNRUN.Nanog.5p.counts.bigWig
    ln -s $basedir/H1_Nanog-GSM3677838/H1_Nanog-GSM3677838.2x25mers.hg38-no_haps.unique.maxFL120bp.5p.counts.bigWig CUTNRUN.Nanog.maxFL120bp.5p.counts.bigWig
    ln -s $basedir/H1_Nanog-GSM3677838/H1_Nanog-GSM3677838.2x25mers.hg38-no_haps.unique.minFL150bp.5p.counts.bigWig CUTNRUN.Nanog.minFL150bp.5p.counts.bigWig
    #Compress peak file
    cat $basedir/H1_Nanog-GSM3677838/H1_Nanog-GSM3677838.2x25mers.hg38-no_haps.unique.MACS2.IDR-0.05 | gzip -c > CUTNRUN.Nanog.narrowpeak.gz
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
    ln -s $basedir/ENCODE-H1-ChIP-seq-datasets/H1-POU5F1-ENCSR000BMU.optimal_idr_thresholded.bed.gz ChIPseq.POU5F1.idroptimal.narrowpeak.gz
    #for NANOG
    ln -s $basedir/ENCODE-H1-ChIP-seq-datasets/H1-NANOG-ENCSR000BMT-ENCFF309DIF.5p.counts.minus.bigWig ChIPseq.NANOG.ENCFF309DIF.5p.counts.minus.bigWig
    ln -s $basedir/ENCODE-H1-ChIP-seq-datasets/H1-NANOG-ENCSR000BMT-ENCFF309DIF.5p.counts.plus.bigWig ChIPseq.NANOG.ENCFF309DIF.5p.counts.plus.bigWig
    ln -s $basedir/ENCODE-H1-ChIP-seq-datasets/H1-NANOG-ENCSR000BMT-ENCFF805SCJ.5p.counts.minus.bigWig ChIPseq.NANOG.ENCFF805SCJ.5p.counts.minus.bigWig
    ln -s $basedir/ENCODE-H1-ChIP-seq-datasets/H1-NANOG-ENCSR000BMT-ENCFF805SCJ.5p.counts.plus.bigWig ChIPseq.NANOG.ENCFF805SCJ.5p.counts.plus.bigWig
    #IDR optimal file for the Nanog ENCSR000BMT experiment
    #For NANOG, the experiment was put through the ENCODE pipeline,
    # so we will just download the IDR optimal file from there
    wget https://www.encodeproject.org/files/ENCFF794GVQ/@@download/ENCFF794GVQ.bed.gz -O ChIPseq.NANOG.idroptimal.narrowpeak.gz
    #For the ChIP-seq control, we will download the .bam file associated with the control
    # for H1ESCs listed on the ENCODE experiment pages for the Nanog and Pou5f1
    # experiments. We will convert the .bam file to a bigwig file.
    wget https://www.encodeproject.org/files/ENCFF948LKH/@@download/ENCFF948LKH.bam -O ChIPseq.control.bam  
}

gather_relevant_files /data/chipseq_vs_newtech/2019-06-26-CUTnRUN-EChO



#!/usr/bin/env nextflow

/*  Nextflow pipeline for create individual phylogenetic trees for each lineage in a WGS dataset. The 
*   files required for Snippy are collected and correctly formats the files required for Snippy into 
*   folder for each separate WGS cluster or lineage.  Snippy is then run on each group separately, 
*   followed by iqtree to generate phyolgeny. 
*
*   written by ellisrichardj
*
*   Version 0.0.1   01/07/20    Initial version - getting started
*
*/

params.ClusterCsv = "$PWD/*_AssignedWGSCluster_*.csv"
params.RefFile = "$HOME/MyScripts/BovTB-nf/references/Mycbovis-2122-97_LT708304.fas"
params.Outdir = 'GroupData'

RefFile = file(params.RefFile)

// This section creates a tuple of samples assigned to each cluster
Channel
    .fromPath( params.ClusterCsv )
    .splitCsv(header: true)
    .filter { row -> row.group =~ /^B.*/ }
    .map { row -> tuple ( row.group, row.Sample ) }
    .into { ClusterSplit1 ; ClusterSplit2 }

//ClusterSplit2.println()

// The groupSamples process collects the required files for each sample and assigns them to a group
process groupSamples {
    errorStrategy 'finish'
    tag "$Sample"

    publishDir 'GroupData', mode: 'copy'

    input:
        set group, val(Sample) from ClusterSplit1

    output:
        set group, file ("${group}/${Sample}/") into TreeGroup 
    
    """
    mkdir -p ${group}/${Sample}
    sed '/^>/ s/.*/>LT708304-Mycobacteriumbovis-AF2122-97/' /$PWD/consensus/${Sample}_consensus.fas > ${group}/${Sample}/${Sample}--${group}_consensus.aligned.fa
    cat /$PWD/snpTables/${Sample}_snps.tab > ${group}/${Sample}/${Sample}--${group}_snps.tab
    gunzip -c /$PWD/vcf/${Sample}.norm.vcf.gz > ${group}/${Sample}/${Sample}--${group}_consensus.vcf
    """
}

//This makes a tuple for each group (cluster), listing path to data for each sample assigned to a particular group
TreeGroup
    .groupTuple()
    .set{ eachTree }

//The PlantTrees process uses Snippy to generate core alignments for all samples in each group
process PlantTrees {
    tag "$group"
    publishDir 'GroupData', mode: 'copy', saveAs: { filename -> "$group/$filename" }

    input: 
    set group, file (samplePaths) from eachTree

    output:
    set group, file ("${group}_core.txt"), file ("${group}_core.aln") into alignments
    
    """
    ~/Tools/snippy-4.6.0/bin/snippy-core --ref ${RefFile} $samplePaths
    mv core.aln ${group}_core.aln
    mv core.txt ${group}_core.txt
    """
}

// The GrowTrees process runs standard iqtree command on the alignment for each of the groups
process GrowTrees {
    tag "$group"
    errorStrategy 'ignore'

    publishDir 'GroupData', mode: 'copy', saveAs: { filename -> "$group/$filename" }

    input:
    set group, file ("${group}_core.txt"), file ("${group}_core.aln") from alignments

    output:
    set group, file ("${group}_core.aln.treefile") into Trees

    """
    ~/Tools/iqtree-2.0.4-Linux/bin/iqtree2 -s ${group}_core.aln
    """
}

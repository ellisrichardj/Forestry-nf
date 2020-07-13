#!/usr/bin/env nextflow

/*  Nextflow pipeline for create individual phylogenetic trees for each lineage in a WGS dataset
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

Channel
    .fromPath( params.ClusterCsv )
    .splitCsv(header: true)
    .map { row -> tuple ( row.group, row.Sample ) }
    .into { ClusterSplit1 ; ClusterSplit2 }


//  The groupSamples process collects and correctly formats the files required for Snippy into a 
//  folder for each separate WGS cluster or lineage.  Snippy is then run on each group separately,
//  followed by iqtree to generate phyolgeny. 

//The groupSamples process collects the required files for each sample and assigns them to a group
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

//The GrowTrees process uses Snippy to generate core alignments for all samples in each group
process GrowTrees {
    tag "$group"
    publishDir 'GroupData', mode: 'copy', saveAs: { filename -> "$group/$filename" }

    input: 
    set group, file (samplePaths) from eachTree

    output:
    set group, file ('*_core.txt'), file ('*_core.aln') into alignments
    
    """
    ~/Tools/snippy-4.6.0/bin/snippy-core --ref ${RefFile} $samplePaths
    mv core.aln ${group}_core.aln
    mv core.txt ${group}_core.txt
    """

}

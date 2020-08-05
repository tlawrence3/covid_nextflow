#!/usr/bin/env nextflow
params.reference="$launchDir/NC_045512.fa"
params.genomes="$launchDir/gisaid_hcov-19_2020_04_28_23.fasta"

reference = file(params.reference)
genomes = file(params.genomes)
mat_peptides = Channel.of( ["nsp1", "266..805"], ["nsp2", "806..2719"], 
                           ["nsp3", "2720..8554"], ["nsp4", "8555..10054"],
                           ["nsp5", "10055..10972"], ["nsp6", "10973..11842"],
                           ["nsp7", "11843..12091"], ["nsp8", "12092..12685"],
                           ["nsp9", "12686..13024"], ["nsp10", "13025..13441"],
                           ["nsp12", "13442..13468,13468..16236"], 
                           ["nsp13", "16237..18039"], ["nsp14", "18040..19620"],
                           ["nsp15", "19621..20658"], ["nsp16", "20659..21552"],
                           ["nsp11", "13442..13480"], ["spike", "21563..25384"],
                           ["orf3a", "25393..26220"], ["orf4", "26245..26472"],
                           ["orf5", "26523..27191"], ["orf6", "27202..27387"],
                           ["orf7a", "27394..27759"], ["orf7b", "27756..27887"],
                           ["orf8", "27894..28259"], ["orf9", "28274..29533"] )

mat_peptides.into { mat_peptides_AA; mat_peptides_cds }

process alignReference {
  publishDir "$launchDir/results/genome_alignment", mode: 'copy'
  
  output:
  path "${genomes.baseName}.aln" into reference_aligned

  shell:
  """
     mafft --auto --thread -1 --keeplength --addfragments !{genomes} !{reference} > !{genomes.baseName}.aln
  """
}

process filterAlignment {
    conda '/home/cades/anaconda3/envs/cov_align' 
    publishDir "$launchDir/results/genome_alignment", mode: 'copy'
    
    input:
    path align_file from reference_aligned

    output:
    file "${align_file.baseName}.filtered.aln" into filtered_align 
    
    script:
    """
    python ${launchDir}/filtering.py ${align_file} ${align_file.baseName}.filtered.aln
    """
}

process extract_genes {
    
    publishDir "$launchDir/results/AA_sequences/initial", mode: 'copy'

    input:
    tuple gene_name, coords from mat_peptides_AA
    each file(genomes) from filtered_align

    output:
    file "${gene_name}.${genomes.baseName}.faa" into proteins
    
    script:
    """
    fascut ${coords} ${genomes} | fastr --degap | fasxl | fassub '\\-xl0' '' > ${gene_name}.${genomes.baseName}.faa
    """
}

process extract_cds {
    publishDir "$launchDir/results/cds_sequences/initial", mode: 'copy'
    
    input:
    tuple gene_name, coords from mat_peptides_cds
    each file(genomes) from filtered_align
    
    output:
    file "${gene_name}.${genomes.baseName}.cds" into cds

    script:
    """
    fascut ${coords} ${genomes} | fastr --degap > ${gene_name}.${genomes.baseName}.cds
    """ 
}

process filter_protein {
    
    conda '/home/cades/anaconda3/envs/cov_align'
    publishDir "$launchDir/results/AA_sequences/filtered", mode: 'copy'
    
    input:
    file protein from proteins
    
    output:
    file "${protein.baseName}.final.faa" into filtered_proteins

    script:
    if( "${protein}" =~ /nsp/ )
        """
        python ${launchDir}/filteringAA.py ${protein} ${protein.baseName}.final.faa 0
        """
    
    else if( "${protein}" =~ /spike/ )
        """
        python ${launchDir}/filteringAA.py ${protein} ${protein.baseName}.final.faa 1
        """
    
    else( "${protein}" =~ /orf/ )
        """
        python ${launchDir}/filteringAA.py ${protein} ${protein.baseName}.final.faa 1
        """
}

process align_protein {
    
    publishDir "$launchDir/results/AA_sequences/alignment", mode: 'copy'

    input:
    file protein from filtered_proteins

    output:
    file "${protein.baseName}.final.aligned.faa" into aligned_proteins    

    script:
    """
    mafft --auto --thread 4 ${protein} > ${protein.baseName}.final.aligned.faa  
    """    
}

process calc_entropy {
    
    conda '/home/cades/anaconda3/envs/cov_align'
    publishDir "$launchDir/results", mode: 'copy'

    input:
    file('*') from aligned_proteins.collect()

    output:
    file "entropy.csv" into entropy_results

    script:
    """
    python ${launchDir}/entropy.py * > entropy.csv
    """
}

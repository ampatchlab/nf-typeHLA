/*
vim: syntax=groovy
-*- mode: groovy;-*-
 *
 * ampatchlab/nf-typeHLA: HLA typing Nextflow pipeline
 *
 * Copyright (C) 2020 QIMR Berghofer Medical Research Institute
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/*
 * Params
 */

params.publish_dir = './results'
params.publish_everything = false
params.publish_mode = 'copy'

params.publish_bwamem = false
params.publish_hla_typing = false


/*
 * Processes
 */

process bwamem {

    tag { sample == readgroup ? sample : "${sample}:${readgroup}" }

    label 'bwakit'

    publishDir(
        path: "${params.publish_dir}/bwamem/${sample}",
        enabled: params.publish_everything || params.publish_bwamem,
        mode: params.publish_mode,
    )

    input:
    tuple sample, readgroup, path(reads)
    path bwa_index

    output:
    tuple sample, path("${readgroup}.aln.bam"), emit: alignments
    tuple sample, path("${readgroup}.hla.HLA-*.fq"), emit: hla_reads

    script:
    def task_cpus = task.cpus > 1 ? task.cpus - 1 : task.cpus

    def idxbase = bwa_index.first().baseName
    def fastq_files = reads.collect { /"$it"/ }.join(' ')

    """
    bwa mem \\
        -t ${task_cpus} \\
        -R '@RG\\tID:${readgroup}\\tSM:${sample}' \\
        "${idxbase}" \\
        ${fastq_files} |
    bwa-postalt.js \\
        -p "${readgroup}.hla" \\
        "${idxbase}.alt" |
    samtools view \\
        -1 \\
        -o "${readgroup}.aln.bam" \\
        -
    """
}

process hla_typing {

    tag { sample }

    label 'bwakit'

    publishDir(
        path: "${params.publish_dir}/hla_typing/${sample}",
        enabled: params.publish_everything || params.publish_hla_typing,
        mode: params.publish_mode,
    )

    input:
    tuple sample, path("staged/fastq/*")
    path("staged/resource_bundle")

    output:
    tuple sample, path("*/*.{fq,mag,sam,tsv}")

    shell:
    '''
    mkdir -p resource_bundle/HLA-ALT-idx
    pushd resource_bundle
    ln -s ../staged/resource_bundle/HLA-ALT-exons.bed
    ln -s ../staged/resource_bundle/HLA-CDS.fa
    pushd HLA-ALT-idx
    ln -s ../../staged/resource_bundle/HLA-ALT-idx/*.fa.gz .
    popd
    popd
    for fasta in resource_bundle/HLA-ALT-idx/*.fa.gz ; do
        bwa index "${fasta}"
    done

    declare -A genes

    for fq in staged/fastq/*.fq ; do
        [[ "${fq}" =~ \\.(HLA-[A-Z]+[0-9]*)\\.fq$ ]] || continue
        genes["${BASH_REMATCH[1]}"]=1
    done

    for gene in "${!genes[@]}" ; do
        echo >&2 "*** Processing gene ${gene}..."

        prefix="!{sample}.${gene}"

        mkdir -p "${prefix}/tmp"
        pushd "${prefix}"

        cat ../staged/fastq/*.${gene}.fq > "${prefix}.fq"
        ln -s ../resource_bundle/* .

        touch "${prefix}.mag"
        touch "${prefix}.sam"
        touch "${prefix}.tsv"

        echo >&2 "** De novo assembling..."

        readlength="$(seqtk comp "${prefix}.fq" | awk '{ s+=$2 } END { printf("%.0f", NR ? s/NR : 0) }')"
        [ "${readlength}" -ge 70 ] || { popd; continue; }

        fermi2.pl unitig -p "./tmp/${prefix}.tmp" -l "${readlength}" -t 1 "${prefix}.fq" | make -f-

        gunzip -dc "./tmp/${prefix}.tmp.mag.gz" > "${prefix}.mag"
        [ -s "${prefix}.mag" ] || { popd; continue; }

        echo >&2 "** Selecting contigs overlapping target exons..."

        for idxbase in HLA-ALT-idx/*.fa.gz ; do
            bwa mem -B1 -O1 -E1 "${idxbase}" "${prefix}.mag" |
            samtools sort -o "./tmp/${prefix}.tmp.ALT.$(basename "${idxbase}" '.fa.gz').bam" -
        done

        samtools merge -O sam "./tmp/${prefix}.tmp.ALT.sam" ./tmp/*.bam

        typeHLA-selctg.js "${gene}" HLA-ALT-exons.bed "./tmp/${prefix}.tmp.ALT.sam" |
        seqtk subseq "${prefix}.mag" - > "./tmp/${prefix}.tmp.fq"

        echo >&2 "** Mapping exons to de novo contigs..."

        bwa index -p "./tmp/${prefix}.tmp" "./tmp/${prefix}.tmp.fq"
        seqtk comp HLA-CDS.fa | cut -f1 | grep "^${gene}" | seqtk subseq HLA-CDS.fa - |
        bwa mem -aD.1 "./tmp/${prefix}.tmp" - > "${prefix}.sam"

        echo >&2 "** Typing..."

        typeHLA.js "${prefix}.sam" > "${prefix}.gt"
        grep "^GT" "${prefix}.gt" | sed "s/^GT/!{sample}/" > "${prefix}.tsv"

        popd
    done
    '''
}

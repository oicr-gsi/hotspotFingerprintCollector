# hotspotFingerprintCollector

fingerprintsCollector, workflow that generates aligns reads to reference, then creates fingerprints using variousmethods. Output are fingerprints of various types and coverage statistics from the alignment
##

## Overview

## Dependencies

* [gatk 4.2.0.0](https://gatk.broadinstitute.org)
* [tabix 0.2.6](http://www.htslib.org)



## Usage

### Cromwell
```
java -jar cromwell.jar run hotspotFingerprintCollector.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputType`|String|one of either fastq or bam
`aligner`|String|aligner to use for fastq input, either bwa or star
`markDups`|String|should the alignment be duplicate marked?, generally yes
`outputFileNamePrefix`|String|Optional output prefix for the output
`refFasta`|String|Path to the reference fasta file
`hotspots`|String|hotspot file for generateFingerprints
`downsample.modules`|String|Names and versions of modules
`bwaMem.runBwaMem_bwaRef`|String|The reference genome to align the sample with by BWA
`bwaMem.runBwaMem_modules`|String|Required environment modules
`markDuplicates.modules`|String|Names and versions of modules
`assessCoverage.modules`|String|Names and versions of modules
`generateFingerprint.modules`|String|Names and versions of modules
`runFinCreator.modules`|String|Names and versions of modules


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`fastqR1`|File?|None|fastq file for read 1
`fastqR2`|File?|None|fastq file for read 2
`bam`|File?|None|bam file
`bamIndex`|File?|None|bam index file
`maxReads`|Int|0|The maximum number of reads to process; if set, this will sample the requested number of reads


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`downsample.jobMemory`|Int|8|memory allocated for Job
`downsample.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`bwaMem.adapterTrimmingLog_timeout`|Int|48|Hours before task timeout
`bwaMem.adapterTrimmingLog_jobMemory`|Int|12|Memory allocated indexing job
`bwaMem.indexBam_timeout`|Int|48|Hours before task timeout
`bwaMem.indexBam_modules`|String|"samtools/1.9"|Modules for running indexing job
`bwaMem.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`bwaMem.bamMerge_timeout`|Int|72|Hours before task timeout
`bwaMem.bamMerge_modules`|String|"samtools/1.9"|Required environment modules
`bwaMem.bamMerge_jobMemory`|Int|32|Memory allocated indexing job
`bwaMem.runBwaMem_timeout`|Int|96|Hours before task timeout
`bwaMem.runBwaMem_jobMemory`|Int|32|Memory allocated for this job
`bwaMem.runBwaMem_threads`|Int|8|Requested CPU threads
`bwaMem.runBwaMem_addParam`|String?|None|Additional BWA parameters
`bwaMem.adapterTrimming_timeout`|Int|48|Hours before task timeout
`bwaMem.adapterTrimming_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.adapterTrimming_addParam`|String?|None|Additional cutadapt parameters
`bwaMem.adapterTrimming_modules`|String|"cutadapt/1.8.3"|Required environment modules
`bwaMem.slicerR2_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR2_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR2_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.slicerR1_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR1_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR1_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.countChunkSize_timeout`|Int|48|Hours before task timeout
`bwaMem.countChunkSize_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.numChunk`|Int|1|number of chunks to split fastq file [1, no splitting]
`bwaMem.trimMinLength`|Int|1|minimum length of reads to keep [1]
`bwaMem.trimMinQuality`|Int|0|minimum quality of read ends to keep [0]
`bwaMem.adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|adapter sequence to trim from read 1 [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]
`bwaMem.adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|adapter sequence to trim from read 2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]
`star.indexBam_timeout`|Int|48|hours before task timeout
`star.indexBam_modules`|String|"picard/2.19.2"|modules for running indexing job
`star.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`star.runStar_timeout`|Int|72|hours before task timeout
`star.runStar_jobMemory`|Int|64|Memory allocated for this job
`star.runStar_threads`|Int|6|Requested CPU threads
`star.runStar_peOvMMp`|Float|0.1|maximum proportion of mismatched bases in the overlap area
`star.runStar_chimSegmentReadGapMax`|Int|3|maximum gap in the read sequence between chimeric segments
`star.runStar_peOvNbasesMin`|Int|10|minimum number of overlap bases to trigger mates merging and realignment
`star.runStar_chimOutJunForm`|Int?|None|flag to add metadata to chimeric junction output for functionality with starFusion - 1 for metadata, 0 for no metadata
`star.runStar_chimNonchimScoDMin`|Int|10|to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value
`star.runStar_chimMulmapNmax`|Int|50|maximum number of chimeric multi-alignments
`star.runStar_chimScoreSeparation`|Int|1|minimum difference (separation) between the best chimeric score and the next one
`star.runStar_chimScoJunNonGTAG`|Int|0|penalty for a non-GTAG chimeric junction
`star.runStar_chimMulmapScoRan`|Int|3|the score range for multi-mapping chimeras below the best chimeric score
`star.runStar_alignIntMax`|Int|100000|maximum intron size
`star.runStar_alignMatGapMax`|Int|100000|maximum gap between two mates
`star.runStar_alignSJDBOvMin`|Int|10|minimum overhang for annotated spliced alignments
`star.runStar_chimJunOvMin`|Int|10|minimum overhang for a chimeric junction
`star.runStar_chimSegmin`|Int|10|minimum length of chimeric segment length
`star.runStar_multiMax`|Int|-1|multiMax parameter for STAR
`star.runStar_saSparsed`|Int|2|saSparsed parameter for STAR
`star.runStar_uniqMAPQ`|Int|255|Score for unique mappers
`star.runStar_chimScoreDropMax`|Int|30|max drop (difference) of chimeric score (the sum of scores of allchimeric segments) from the read length
`star.runStar_outFilterMultimapNmax`|Int|50|max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped
`star.runStar_modules`|String|"star/2.7.6a hg38-star-index100/2.7.6a"|modules for running STAR
`star.runStar_addParam`|String?|None|Additional STAR parameters
`star.runStar_genereadSuffix`|String|"ReadsPerGene.out"|ReadsPerGene file suffix
`star.runStar_chimericjunctionSuffix`|String|"Chimeric.out"|Suffix for chimeric junction file
`star.runStar_transcriptomeSuffix`|String|"Aligned.toTranscriptome.out"|Suffix for transcriptome-aligned file
`star.runStar_starSuffix`|String|"Aligned.sortedByCoord.out"|Suffix for sorted file
`star.runStar_genomeIndexDir`|String|"$HG38_STAR_INDEX100_ROOT/"|Path to STAR index
`markDuplicates.jobMemory`|Int|24|memory allocated for Job
`markDuplicates.overhead`|Int|6|memory allocated to overhead of the job other than used in markDuplicates command
`markDuplicates.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`assessCoverage.jobMemory`|Int|8|memory allocated for Job
`assessCoverage.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`generateFingerprint.jobMemory`|Int|8|memory allocated for Job
`generateFingerprint.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`generateFingerprint.stdCC`|Float|30.0|tunable standard call confidence score parameter for gatk, default is 30
`runFinCreator.chroms`|Array[String]|["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]|Canonical chromosomes in desired order (used for soting lines in .fin file)
`runFinCreator.timeout`|Int|10|Timeout in hours, needed to override imposed limits
`runFinCreator.jobMemory`|Int|8|memory allocated for Job


### Outputs

Output | Type | Description | Labels
---|---|---|---
`coverage`|File|output from samtools coverage, with per chromosome metrics|vidarr_label: coverage 
`json`|File|metrics in json format, currently only the mean coverage for the alignment|vidarr_label: json 
`hotspotVcf`|File|vcf file from GATK haplotype caller on provided hotspots|vidarr_label: hotspotVcf 
`hotspotTbi`|File|index for the vcf hotspot fingerprint|vidarr_label: hotspotTbi 
`hotspotFin`|File|fingerprint Fin file for the provided hotspots|vidarr_label: hotspotFin 


## Commands
This section lists command(s) run by hotspotFingerprintsCollector workflow
 
* Running hotspotFingerprintsCollector
 

 
### Commands for downsampling
 
```
  set -euo pipefail
  
  seqtk sample -s 100 INPUT_FASTQ_R1 maxReads > OUTPUT_FATSQ_R1m
  gzip OUTPUT_FATSQ_R1m
  
  seqtk sample -s 100 INPUT_FASTQ_R2 maxReads > OUTPUT_FATSQ_R2m
  gzip OUTPUT_FATSQ_R2m
```

### Command for markDuplicates 

```
   set -euo pipefail
   $GATK_ROOT/bin/gatk --java-options "-Xmx [JOB_MEMORY] - [OVERHEAD]G" MarkDuplicates \
                       -I INPUT_Bam \
                       --METRICS_FILE OUTPUT_FILE_NAME_PREFIX.dupmetrics \
                       --VALIDATION_STRINGENCY SILENT \
                       --CREATE_INDEX true \
                       -O OUTPUT_FILE_NAME_PREFIX.dupmarked.bam
```
 
### Commands for generateFingerprint

```
   set -euo pipefail
   $SAMTOOLS_ROOT/bin/samtools coverage INPUT_BAM > OUTPUT_FILE_NAME_PREFIX.coverage.txt
   cat OUTPUT_FILE_NAME_PREFIX.coverage.txt | grep -P "^chr\d+\t|^chrX\t|^chrY\t" | awk '{ space += ($3-$2)+1; bases += $7*($3-$2);} END { print bases/space }' | awk '{print "{\"mean coverage\":" $1 "}"}' > OUTPUT_FILE_NAME_PREFIX.json
 
 
  set -euo pipefail
  $GATK_ROOT/bin/gatk HaplotypeCaller \
                     -R refFasta \
                     -I INPUT_BAM \
                     -O OUTPUT_FILE_NAME_PREFIX.hotspots.vcf \
                     --read-filter CigarContainsNoNOperator \
                     --stand-call-conf stdCC \
                     -L hotspots
 
  $TABIX_ROOT/bin/bgzip -c OUTPUT_FILE_NAME_PREFIX.hotspots.vcf > OUTPUT_FILE_NAME_PREFIX.hotspots.vcf.gz
  $TABIX_ROOT/bin/tabix -p vcf OUTPUT_FILE_NAME_PREFIX.hotspots.vcf.gz 
  
  $GATK_ROOT/bin/gatk DepthOfCoverage \
                      -R refFasta \
                      -I INPUT_BAM \
                      -O OUTPUT_FILE_NAME_PREFIX \
                      --read-filter CigarContainsNoNOperator \
                      -L hotspots 
  
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_

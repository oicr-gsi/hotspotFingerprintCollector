version 1.0

import "imports/pull_bwaMem.wdl" as bwaMem
import "imports/pull_star.wdl" as star

struct InputGroup {
  File fastqR1
  File fastqR2
  String readGroup
}

workflow fingerprintsCollector {
   input {
        File? fastqR1
        File? fastqR2
        File? bam
        File? bamIndex
        String inputType
        String aligner
        String markDups
        String outputFileNamePrefix
        String refFasta
        String haplotypeMap
        Int maxReads = 0
        File? hotspots
   }
   
   parameter_meta {
        fastqR1: "fastq file for read 1"
        fastqR2: "fastq file for read 2"
        bam: "bam file"
        bamIndex: "bam index file"
        inputType: "one of either fastq or bam"
        aligner : "aligner to use for fastq input, either bwa or star"
        markDups : "should the alignment be duplicate marked?, generally yes"
        outputFileNamePrefix: "Optional output prefix for the output"
        refFasta: "Path to the reference fasta file"
        haplotypeMap: "Path to the gzipped hotspot vcf file"
        maxReads: "The maximum number of reads to process; if set, this will sample the requested number of reads"
		hotspots: "hotspot file for generateFingerprints"
   }

   if(inputType=="fastq" && defined(fastqR1) && defined(fastqR2)){
     if(maxReads>0){
      call downsample{
        input:
          fastqR1 = select_first([fastqR1]),
          fastqR2 = select_first([fastqR2]),
          maxReads = maxReads
      }
     }

     if(aligner=="bwa"){
       call bwaMem.bwaMem {
         input:
           fastqR1 = select_first([downsample.fastqR1mod,fastqR1]),
           fastqR2 = select_first([downsample.fastqR2mod,fastqR2]),
           outputFileNamePrefix = outputFileNamePrefix,
           readGroups = "'@RG\\tID:CROSSCHECK\\tSM:SAMPLE'",
           doTrim = false
        }
      }

      if(aligner=="star"){
       InputGroup starInput = { 
         "fastqR1": select_first([downsample.fastqR1mod,fastqR1]),
         "fastqR2": select_first([downsample.fastqR2mod,fastqR2]),
         "readGroup": "ID:CROSSCHECK SM:SAMPLE"
       }
       call star.star { 
         input:
           inputGroups = [ starInput ],
           outputFileNamePrefix = outputFileNamePrefix,
           runStar_chimOutType = "Junctions"
       }
     }
   }
  if(markDups=="true"){
    call markDuplicates {
      input :
        inputBam = select_first([bwaMem.bwaMemBam,star.starBam,bam]),
        inputBai = select_first([bwaMem.bwaMemIndex,star.starIndex,bamIndex]),
        outputFileNamePrefix = outputFileNamePrefix 
    }
  }  

   call assessCoverage {
     input:
        inputBam = select_first([markDuplicates.bam,bwaMem.bwaMemBam,star.starBam,bam]),
        inputBai = select_first([markDuplicates.bamIndex,bwaMem.bwaMemIndex,star.starIndex,bamIndex]),
        outputFileNamePrefix = outputFileNamePrefix  
   }
   
   call extractFingerprint {
     input:
        inputBam = select_first([markDuplicates.bam,bwaMem.bwaMemBam,star.starBam,bam]),
        inputBai = select_first([markDuplicates.bamIndex,bwaMem.bwaMemIndex,star.starIndex,bamIndex]),
        haplotypeMap = haplotypeMap,
        refFasta = refFasta,
        outputFileNamePrefix = outputFileNamePrefix
    }

   call generateFingerprint {
     input:
        inputBam = select_first([markDuplicates.bam,bwaMem.bwaMemBam,star.starBam,bam]),
        inputBai = select_first([markDuplicates.bamIndex,bwaMem.bwaMemIndex,star.starIndex,bamIndex]),
        hotspots = select_first([hotspots]),
        refFasta = refFasta,
        outputFileNamePrefix = outputFileNamePrefix
   }

   call runFinCreator { 
     input: 
        inputVcf = generateFingerprint.vcf, 
        inputCoverage = generateFingerprint.depth, 
        sampleID = outputFileNamePrefix,
        hotspotSNPs = select_first([hotspots])
   }
   


   output {
      File outputVcf = extractFingerprint.vgz
      File outputTbi = extractFingerprint.tbi
      File coverage = assessCoverage.coverage
      File json = assessCoverage.json
      File hotspotVcf = generateFingerprint.vgz
      File hotspotTbi = generateFingerprint.tbi
      File hotspotFin = runFinCreator.fin
     }

    meta {
     author: "Lawrence Heisler"
     email: "lawrence.heisler@oicr.on.ca"
     description: "crosscheckFingerprintsCollector, workflow that generates genotype fingerprints using gatk ExtractFingprint.  Output are vcf files that can be proccessed through gatk Crosscheck fingerprints\n##"
     dependencies: [
      {
        name: "gatk/4.2.0.0",
        url: "https://gatk.broadinstitute.org"
      },
      {
        name: "tabix/0.2.6",
        url: "http://www.htslib.org"
      }
     ]
     output_meta: {
       outputVcf: "the crosscheck fingerprint, gzipped vcf file",
       outputTbi: "index for the vcf fingerprint",
       coverage : "output from samtools coverage, with per chromosome metrics",
       json : "metrics in json format, currently only the mean coverage for the alignment",
       hotspotVcf : "vcf file from GATK haplotype caller on provided hotspots",
       hotspotTbi : "index for the vcf hotspot fingerprint",
       hotspotFin : "fingerprint Fin file for the provided hotspots"

     }
  }

}
# ==========================================
#  configure and run extractFingerprintsCollector
# ==========================================


task extractFingerprint {
input {
 File inputBam
 File inputBai
 String modules
 String refFasta
 String outputFileNamePrefix
 String haplotypeMap
 Int jobMemory = 8
 Int timeout = 24
}
parameter_meta {
 inputBam: "input .bam file"
 inputBai: "index of the input .bam file"
 refFasta: "Path to reference FASTA file"
 outputFileNamePrefix: "prefix for making names for output files"
 haplotypeMap: "Hotspot SNPs are the locations of variants used for genotyping"
 jobMemory: "memory allocated for Job"
 modules: "Names and versions of modules"
 timeout: "Timeout in hours, needed to override imposed limits"
}

command <<<
  set -euo pipefail

 $GATK_ROOT/bin/gatk ExtractFingerprint \
                    -R ~{refFasta} \
                    -H ~{haplotypeMap} \
                    -I ~{inputBam} \
                    -O ~{outputFileNamePrefix}.vcf
 $TABIX_ROOT/bin/bgzip -c ~{outputFileNamePrefix}.vcf > ~{outputFileNamePrefix}.vcf.gz
 $TABIX_ROOT/bin/tabix -p vcf ~{outputFileNamePrefix}.vcf.gz 


>>>

 runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
 }

 output {
  File vcf = "~{outputFileNamePrefix}.vcf"
  File vgz = "~{outputFileNamePrefix}.vcf.gz"
  File tbi = "~{outputFileNamePrefix}.vcf.gz.tbi"
 }
}



# ==========================================
#  downsample the fastq files to the first N reads
# ==========================================


task downsample {
 input{
  File fastqR1
  File fastqR2
  Int maxReads
  String modules
  Int jobMemory = 8
  Int timeout = 24
 }
 
 parameter_meta {
  fastqR1 : "Read1 fastq file"
  fastqR2 : "Read2 fastq file"
  maxReads : "the maximum number of reads to use"
  jobMemory: "memory allocated for Job"
  modules: "Names and versions of modules"
  timeout: "Timeout in hours, needed to override imposed limits"
 }
 
 String fastqR1m = basename(fastqR1,".fastq.gz") + ".mod.fastq"
 String fastqR2m = basename(fastqR2,".fastq.gz") + ".mod.fastq"
  
command <<<
 set -euo pipefail
 
 seqtk sample -s 100 ~{fastqR1} ~{maxReads} > ~{fastqR1m}
 gzip ~{fastqR1m}
 
 seqtk sample -s 100 ~{fastqR2} ~{maxReads} > ~{fastqR2m}
 gzip ~{fastqR2m}
>>>

 runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
 }

 output {
  File fastqR1mod = "~{fastqR1m}.gz"
  File fastqR2mod = "~{fastqR2m}.gz"
 }    
}


# ==========================================
#  Duplicate Marking
# ==========================================



task markDuplicates {
 input{
  File inputBam
  File inputBai
  String modules
  String outputFileNamePrefix
  Int jobMemory = 8
  Int timeout = 24
 }
 parameter_meta {
  inputBam: "input .bam file"
  inputBai: "index of the input .bam file"
  outputFileNamePrefix: "prefix for making names for output files"  
  jobMemory: "memory allocated for Job"
  modules: "Names and versions of modules"
  timeout: "Timeout in hours, needed to override imposed limits"
 }
 
command <<<
  set -euo pipefail
  $GATK_ROOT/bin/gatk MarkDuplicates \
                      -I ~{inputBam} \
                      --METRICS_FILE ~{outputFileNamePrefix}.dupmetrics \
                      --VALIDATION_STRINGENCY SILENT \
                      --CREATE_INDEX true \
                      -O ~{outputFileNamePrefix}.dupmarked.bam
>>>

 runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
 }

 output {
  File bam = "~{outputFileNamePrefix}.dupmarked.bam"
  File bamIndex = "~{outputFileNamePrefix}.dupmarked.bai"
 }
}

# ==========================================
#  coverage metrics from the bam file used for fingerprint analysis
# ==========================================


 task assessCoverage{
   input{
    File inputBam
    File inputBai
    String modules
    String outputFileNamePrefix
    Int jobMemory = 8
    Int timeout = 24
   }
   parameter_meta {
    inputBam: "input .bam file"
    inputBai: "index of the input .bam file"
    outputFileNamePrefix: "prefix for making names for output files"  
    jobMemory: "memory allocated for Job"
    modules: "Names and versions of modules"
    timeout: "Timeout in hours, needed to override imposed limits"
   } 

command <<<
  set -euo pipefail
  $SAMTOOLS_ROOT/bin/samtools coverage ~{inputBam} > ~{outputFileNamePrefix}.coverage.txt
  cat ~{outputFileNamePrefix}.coverage.txt | grep -P "^chr\d+\t|^chrX\t|^chrY\t" | awk '{ space += ($3-$2)+1; bases += $7*($3-$2);} END { print bases/space }' | awk '{print "{\"mean coverage\":" $1 "}"}' > ~{outputFileNamePrefix}.json
>>>

  runtime {
   memory:  "~{jobMemory} GB"
   modules: "~{modules}"
   timeout: "~{timeout}"
  }

  output {
    File coverage = "~{outputFileNamePrefix}.coverage.txt"
    File json = "~{outputFileNamePrefix}.json"
  }
}



# ==========================================
#  generateFingerprint with gatk haplotype caller
# ==========================================

task generateFingerprint {
input {
 File inputBam
 File inputBai
 String modules
 String refFasta
 String outputFileNamePrefix
 String hotspots
 Int jobMemory = 8
 Int timeout = 24
 Float stdCC = 30.0
}
parameter_meta {
 inputBam: "input .bam file"
 inputBai: "index of the input .bam file"
 refFasta: "Path to reference FASTA file"
 outputFileNamePrefix: "prefix for making names for output files"
 hotspots: "Hotspot SNPs are the locations of variants used for genotyping"
 jobMemory: "memory allocated for Job"
 modules: "Names and versions of modules"
 timeout: "Timeout in hours, needed to override imposed limits"
 stdCC: "tunable standard call confidence score parameter for gatk, default is 30"
}

command <<<
 set -euo pipefail
 $GATK_ROOT/bin/gatk HaplotypeCaller \
                    -R ~{refFasta} \
                    -I ~{inputBam} \
                    -O ~{outputFileNamePrefix}.hotspots.vcf \
                    --read-filter CigarContainsNoNOperator \
                   --stand-call-conf ~{stdCC} \
                    -L ~{hotspots}

 $TABIX_ROOT/bin/bgzip -c ~{outputFileNamePrefix}.hotspots.vcf > ~{outputFileNamePrefix}.hotspots.vcf.gz
 $TABIX_ROOT/bin/tabix -p vcf ~{outputFileNamePrefix}.hotspots.vcf.gz 
 
 $GATK_ROOT/bin/gatk DepthOfCoverage \
                     -R ~{refFasta} \
                     -I ~{inputBam} \
                     -O ~{outputFileNamePrefix} \
                     --read-filter CigarContainsNoNOperator \
                     -L ~{hotspots} 
 
>>>

 runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
 }

 output {
  File vcf = "~{outputFileNamePrefix}.hotspots.vcf"
  File vgz = "~{outputFileNamePrefix}.hotspots.vcf.gz"
  File tbi = "~{outputFileNamePrefix}.hotspots.vcf.gz.tbi"
  File depth = "~{outputFileNamePrefix}.sample_interval_summary"
 }
}



task runFinCreator {
 input {
  File inputVcf
  File inputCoverage
  File hotspotSNPs
  String sampleID
  String modules
  Array[String] chroms = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]
  Int timeout = 10
  Int jobMemory = 8
}

parameter_meta {
 inputVcf: "Input .bam file for analysis sample"
 inputCoverage: "Optional input .bam file for control sample"
 hotspotSNPs: "Hotspot vcf, used as a reference"
 sampleID: "This is used as a prefix for output files"
 jobMemory: "memory allocated for Job"
 modules: "Names and versions of modules"
 chroms: "Canonical chromosomes in desired order (used for soting lines in .fin file)"
 timeout: "Timeout in hours, needed to override imposed limits"
}

command <<<
 python3 <<CODE
 import gzip
 import re

 flags = {"zerocov": "N",  # No coverage
          "nosnp": "M"}  # no SNP, base matches the reference

 # Read lines from the file with reference SNPs
 refFile = gzip.open("~{hotspotSNPs}", mode='r')
 refLines = refFile.readlines()
 refFile.close()

 # Retrieve reference SNPsb"abcde".decode("utf-8")
 refs = {}
 for rLine in refLines:
     refLine = rLine.decode('utf-8')
     if refLine.strip().startswith('#'):
         continue
     temp = refLine.split("\t")
     if not refs.get(temp[0]):
         refs[temp[0]] = {}
     refs[temp[0]][temp[1]] = dict(flag=flags['zerocov'], dbsnp=temp[2], ref=temp[3])

 # Read lines from the file with Depth data
 depthFile = open("~{inputCoverage}", mode='r')
 depthLines = depthFile.readlines()
 depthFile.close()

 # Retrieve Depth values
 for dLine in depthLines:
     if dLine.startswith("Target"):
         continue
     #temp = dLine.split("\t")
     temp = dLine.split(",")
     #coords = temp[0].split(":")
     coords = temp[0].replace("-",":").split(":")
     if len(coords) != 2:
         continue
     if refs[coords[0]].get(coords[1]) and int(temp[1]) > 0:
         refs[coords[0]][coords[1]]['flag'] = flags['nosnp']
     else:
         refs[coords[0]][coords[1]]['flag'] = flags['zerocov']

 # Read lines from the file with Depth data
 callFile = open("~{inputVcf}", mode='r')
 callLines = callFile.readlines()
 callFile.close()

 for cLine in callLines:
     if cLine.startswith('#'):
         continue
     temp = cLine.split("\t")
     if refs[temp[0]] and refs[temp[0]].get(temp[1]):
         variant = temp[4].upper()
         if re.fullmatch('[ACGT]{1}', variant) and temp[4] != temp[3]:
             refs[temp[0]][temp[1]]["flag"] = variant
         else:
             refs[temp[0]][temp[1]]["flag"] = flags['nosnp']

 # Prepare lines for printing
 finLines = []
 for chr in ["~{sep='\",\"' chroms}"]:
     if not refs.get(chr):
         continue
     for start in sorted(refs[chr].keys()):
         snp = [refs[chr][start]['ref']]
         if refs[chr][start]['flag'] == 'M':
             snp.append(refs[chr][start]['ref'])
         else:
             snp.append(refs[chr][start]['flag'])
         if refs[chr][start]['flag'] == 'N':
             snp = [""]
         snpflag = "".join(snp)
         finLines.append("\t".join([chr, start, refs[chr][start]['dbsnp'], snpflag, refs[chr][start]['flag']]))

 # Print into a .fin file
 finHeader = ["CHROM", "POS", "ID", "SNP", "FLAG"]
 f = open("~{sampleID}.fin", "w+")
 f.write('\t'.join(finHeader) + '\n')
 f.write('\n'.join(finLines) + '\n')
 f.close()
 
 CODE
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

 output {
  File fin = "~{sampleID}.fin"
 }
}


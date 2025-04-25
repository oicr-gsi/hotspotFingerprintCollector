# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2025-04-23
- [GRD-795](https://jira.oicr.on.ca/browse/GRD-795) - Expanded built-in documentation (metadata changes only).

## [Unreleased] - 2024-06-25
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - add vidarr labels to outputs (changes to medata only)

## [1.0.0] 2022-03-09
### Added
- Accept fastq files as inputs, also allow bam file as input. Default input type is fastq.
- Generated README.md

### Changed
- Changed repo and workflow name. From "fingerprintCollectror" to "hotspotFingerprintCollector"

### Removed
- The "extractFingerprint" task.

# [1.0.2] - 2022-07-26
- Adding java option Xmx to make sure markDupilcates command running has enough memory 

## [1.0.1] 2022-03-13
### Added 
- None

### Changed
- Moved python code from command line of wdl file to separate python file

### Removed 
- None

## [0.0.9] - 2021-12-16
### Added
- Initiated repo
- Added WDL, subworkflows and imports from crosscheckFingerprintsCollector; will extend to process alignments for other fingerprinting methods
- Added Vidarr files


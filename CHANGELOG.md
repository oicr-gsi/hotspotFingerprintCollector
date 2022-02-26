# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2021-12-16
### Added
- Initiated repo
- Added WDL, subworkflows and imports from crosscheckFingerprintsCollector; will extend to process alignments for other fingerprinting methods
- Added Vidarr files

## [1.1.0] 2022-02-25
### Added
- Accept fastq files as inputs, also allow bam file as input. Default input type is fastq.
- Generated README.md

### Changed
- Changed repo and workflow name. From "fingerprintCollectror" to "hotspotFingerprintCollector"

### Removed
- The "extractFingerprint" task.

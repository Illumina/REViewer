# Repeat Expansion Viewer (REViewer)

REViewer is a tool for visualizing alignments of reads in regions containing
tandem repeats. REViewer requires a BAMlet with graph-realigned reads generated
by [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) and the
corresponding variant catalog.

## License

REViewer is provided under the terms and conditions of the [GPLv3 license](LICENSE.txt).
It relies on several third party packages provided under other open source licenses,
please see [COPYRIGHT.txt](COPYRIGHT.txt) for additional details.

## Installation

The simplest way of obtaining REViewer is by downloading a Linux binary
corresponding to the latest release from the
[Releases page](https://github.com/Illumina/REViewer/releases). The link to the
binary is located in the *Assets* section.

REViewer can also be built from source with CMake.

```shell script
cd REViewer
mkdir build; cd build
cmake ..; make
```

## Usage

REViewer requires output files generated by [ExpansionHunter](https://github.com/Illumina/ExpansionHunter)
v3.0.0 or above along with the
[matching variant catalog](https://github.com/Illumina/ExpansionHunter/blob/master/docs/04_VariantCatalogFiles.md)
file and reference genome.

```shell script
REViewer \
  --reads <BAMlet generated by ExpansionHunter> \
  --vcf <VCF file generated by ExpansionHunter> \
  --reference <FASTA file with reference genome> \
  --catalog <Variant catalog> \
  --locus <Locus to analyze> \
  --output-prefix <Prefix for the output files>
```

Note that the BAMlet generated by ExpansionHunter (`--reads` parameter) must be sorted and indexed.

## How-to guides

- The best way to learn about how to distinguish correctly and incorrectly genotyped
repeats is by going over [examples of read pileups](docs/examples.md).

## Reference documentation

- [A blog post describing the method](https://www.illumina.com/science/genomics-research/reviewer-visualizing-alignments-short-reads-long-repeat.html)
- [Overview of the method and its limitations](docs/method-overview.md)
- [Description of the quality metrics reported by REViewer](docs/metrics.md)

## Companion tools

- [FlipBook](https://github.com/broadinstitute/flipbook) is an image server for
REViewer developed by [Ben Weisburd](https://github.com/bw2). It provides a
convenient way to inspect large quantities of read pileups.

- [Review BAMs](https://gitlab.com/andreassh/review-bams) is script that allows
running REViewer on regular BAM file (by running ExpansionHunter behind the
scenes). It was developed by [Andreas Halman](https://gitlab.com/andreassh).

## Citation

Dolzhenko E, Weisburd B, Garikano K, Rajan Babu IS, and colleagues,
[REViewer: Haplotype-resolved visualization of read alignments in and around tandem repeats](https://www.biorxiv.org/content/10.1101/2021.10.20.465046v1), bioRxiv, 2021

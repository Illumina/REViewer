# Quality metrics

REViewer reports various summary measurements, called **quality metrics**, that
describe key properties of read pileups. Quality metrics make it possible to
automate assessment of large collections of STR genotype calls either (a) by
selecting a series of thresholds to stratify the value of each metric into
"good", "suspicious", and "bad" categories, or (b) by using more flexible
statistical / machine learning approaches.

This document describes the quality metrics reported by REViewer. All quality
metrics are stored in a tab-separated (TSV) file `<output prefix>.metrics.tsv`.
If you have suggestions for additional metrics, please consider [creating an
issue](https://github.com/Illumina/REViewer/issues).

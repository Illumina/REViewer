# Examples of read pileups

Pileups corresponding to correctly genotyped repeats are characterized by a
relatively even read coverage of both alleles (Figure 1, panels 1-3). At the
typical whole-genome sequencing depths (30-60x), each position of a haplotype
sequence is expected to be covered by many reads (15-30), although the coverage
may dip in certain regions due to technical factors like GC bias. For repeats
much shorter than the read length this implies the presence of multiple spanning
reads (Figure 1, both alleles on panel 1 and short allele on panel 2). The
repeats much larger than the read length are expected to contain multiple
in-repeat reads (Figure 1, long allele on panel 1 and both alleles on panel 2).
An expanded allele might not be called correctly if the sequencing depth inside
the repeat is very low compared to the depth of the region surrounding the
repeat (Figure 1, long allele on panel 4). Additionally, the presence of
multiple indels in the alignments of in-repeat reads indicates that the reads
may not be correctly aligned and that the size of the repeat may be
overestimated (Figure 1, panel 5). Finally, a short allele supported by one or
very few spanning reads may not be real. For instance, the short allele depicted
on panel 6 of Figure 1 is supported by just one spanning and one flanking read,
which is less than expected based on the coverage of the surrounding region.
There is also a slight excess of the flanking reads on the long allele of this
repeat. Taken together, these observations suggest that (a) the single spanning
read may be a result of an incorrect alignment and (b) the correct genotype is
likely to be a double expansion. Supplementary Figures S1-6 are the real
examples of situations depicted on panels 1-6 of Figure 1 respectively.

![Examples of read pileups](docs/images/cartoon-examples.png)
**Figure 1: Examples of read pileups.** Pileups corresponding to correctly genotyped
repeats: (1) both repeat alleles are short; (2) one allele is expanded; (3) both
alleles are expanded. Pileups corresponding to incorrectly genotyped repeats:
(4) expanded allele is supported by just one read suggesting that its size is
overestimated; (5) expanded allele is supported by poorly aligning reads (each
containing multiple indels) suggesting that the reads are incorrectly mapped and
that size of the repeat is overestimated; (6) the short allele is supported by
just one spanning read suggesting that this allele is not real and that both
alleles are expanded.



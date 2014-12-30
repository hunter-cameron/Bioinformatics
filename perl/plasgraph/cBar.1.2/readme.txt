cBar: cBar: an algorithm for accurate identification of chromosomal sequence
fragments from metagenome data
Fengfeng Zhou (ffzhou@live.com)


Prerequisition
---

cBar needs the following softwares:

BioPerl       A perl library for bioinformatics computation.
              URL: http://www.bioperl.org/wiki/Main_Page
Weka          A data mining software. cBar has attached Weka 3.6.0
              URL: http://www.cs.waikato.ac.nz/~ml/weka/


Installation
---

Use this command to compile the c codes in your system:

make


Usage
---

You may use this command to process the example genome sequence of Synechococcus elongatus PCC 7942.

> ./cBar.pl Synechococcus_elongatus_PCC_7942.fna Synechococcus_elongatus_PCC_7942.cBar.txt
Formatting the sequences ...  [FormattedSeq:2] [done]
Generating the K-mer profiles ...  [KMer] [Profile] [done]
Predicting ...  [Predicted] [Saved] [Cleared] [done]




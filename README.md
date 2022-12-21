score
=====

A set of tools to handle and convert summary statistics files following the [GWAS-VCF specification](https://github.com/MRCIEU/gwas-vcf-specification). If you use any of these tools in your publication, please cite this website. For any feedback or questions, contact the [author](mailto:giulio.genovese@gmail.com)

We encourage users to adopt the GWAS-VCF specification rather than the [GWAS-SSF specification](https://ebispot.github.io/gwas-blog/new-standard-for-gwas-summary-statistics) promoted by the [GWAS catalog](https://www.ebi.ac.uk/gwas/) as the latter is affected by [issues](https://github.com/EBISPOT/gwas-summary-statistics-standard/issues/4) and furthermore we believe that many common uses are better addressed by using the more general VCF specification. If you are planning to publish your summary statistics, we encourage you to submit them as GWAS-VCF files or as both GWAS-VCF and as GWAS-SSF files. The latter can be generated from the former with the following command
```
(echo -e "chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\teffect_allele_frequency\tp_value";
bcftools query -s SM -f "%CHROM\t%POS\t%ALT\t%REF[\t%ES\t%SE\t%AF\t%LP]\n" gwas-vcf.vcf | \
  sed 's/^chr//;s/^X/23/;s/^Y/24/;s/^MT/25/;s/^M/25/;s/\t\./\tNA/g' | awk -F"\t" -v OFS="\t" '{$8=10^(-$8); print}') > gwas-ssf.tsv
```

<!--ts-->
   * [Usage](#usage)
   * [Installation](#installation)
   * [Column Headers Mappings](#column-headers-mappings)
   * [LDGM-VCF Specification](#ldgm-vcf-specification)
   * [LDGM Matrices](#ldgm-matrices)
   * [Compute polygenic scores](#compute-polygenic-scores)
   * [Convert summary statistics](#convert-summary-statistics)
   * [Liftover VCFs](#liftover-vcfs)
   * [Compute best linear unbiased predictor](#compute-best-linear-unbiased-predictor)
   * [Run meta-analysis](#run-meta-analysis)
   * [Annotation](#annotation)
   * [Plotting](#plotting)
   * [Examples](#examples)
      * [Attention Deficit Hyperactivity Disorder](#attention-deficit-hyperactivity-disorder)
      * [Anxiety Disorder](#anxiety-disorder)
      * [Autism Spectrum Disorder](#autism-spectrum-disorder)
      * [Bipolar Disorder](#bipolar-disorder)
      * [Eating Disorders](#eating-disorders)
      * [Major Depressive Disorder](#major-depressive-disorder)
      * [OCD & Tourette Syndrome](#ocd--tourette-syndrome)
      * [Post Traumatic Stress Disorder](#post-traumatic-stress-disorder)
      * [Schizophrenia](#schizophrenia)
      * [Educational Attainment](#educational-attainment)
      * [Height](#height)
      * [BMI](#bmi)
      * [Smoking](#smoking)
   * [Acknowledgements](#acknowledgements)
<!--te-->

Usage
=====

Polygenic score tool:
```
Usage: bcftools +score [options] <in.vcf.gz> [<score1.gwas.vcf.gz> <score2.gwas.vcf.gz> ...]
Plugin options:
       --use <tag>               FORMAT tag to use to compute allele dosages: GP, AP, HDS, DS, GT, AS
       --summaries <dir|file>    summary statistics files from directory or list from file
       --q-score-thr LIST        comma separated list of p-value thresholds
       --counts                  include SNP counts in the output table
   -o, --output <file.tsv>       write output to a file [standard output]
       --sample-header           output header for sample ID column [SAMPLE]
   -e, --exclude <expr>          exclude sites for which the expression is true
   -f, --apply-filters <list>    require at least one of the listed FILTER strings (e.g. "PASS,.")
   -i, --include <expr>          select sites for which the expression is true
   -r, --regions <region>        restrict to comma-separated list of regions
   -R, --regions-file <file>     restrict to regions listed in a file
       --regions-overlap 0|1|2   Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
   -t, --targets [^]<region>     restrict to comma-separated list of regions. Exclude regions with "^" prefix
   -T, --targets-file [^]<file>  restrict to regions listed in a file. Exclude regions with "^" prefix
       --targets-overlap 0|1|2   Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]
   -s, --samples [^]<list>       comma separated list of samples to include (or exclude with "^" prefix)
   -S, --samples-file [^]<file>  file of samples to include (or exclude with "^" prefix)
       --force-samples           only warn about unknown subset samples

TSV Summary Statistics Options:
   -c, --columns <preset>        column headers from preset (PLINK/PLINK2/REGENIE/SAIGE/METAL/PGS/SSF)
   -C, --columns-file <file>     column headers from tab-delimited file
       --use-variant-id          use variant_id to match variants rather than chromosome and base_pair_location

Examples:
   bcftools +score --use DS -o scores.tsv input.bcf -c PLINK score.assoc
   bcftools +score --use DS -o scores.tsv input.bcf -C colheaders.tsv PGC3_SCZ_wave3_public.clumped.v2.tsv.gz
   bcftools +score --use GT -o scores.tsv --q-score-thr 1e-8,1e-7,1e-6,1e-5,1e-4,0.001,0.01,0.05 input.bcf -c GWAS-SSF PGS000001.txt.gz
   bcftools +score --use DS -o scores.tsv -i 'INFO>0.8 && AF>0.01 && AF<0.99' input.bcf -c GWAS-SSF PGS000001.txt.gz PGS000002.txt.gz
```

Munge summary statistics tool:
```
Usage: bcftools +munge [options] <score.gwas.ssf.tsv>
Plugin options:
   -c, --columns <preset>          column headers from preset (PLINK/PLINK2/REGENIE/SAIGE/METAL/PGS/SSF)
   -C, --columns-file <file>       column headers from tab-delimited file
   -f, --fasta-ref <file>          reference sequence in fasta format
       --fai <file>                reference sequence .fai index
       --set-cache-size <int>      select fasta cache size in bytes
       --iffy-tag <string>         FILTER annotation tag to record whether reference allele could not be determined [IFFY]
   -s, --sample-name <string>      sample name for the phenotype [SAMPLE]
       --ns <float>                number of samples
       --nc <float>                number of cases
       --ne <float>                effective sample size
       --no-version                do not append version and command line to the header
   -o, --output <file>             write output to a file [no output]
   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
       --threads <int>             use multithreading with INT worker threads [0]

Examples:
      bcftools +munge -c PLINK -f human_g1k_v37.fasta -Ob -o score.bcf score.assoc
      bcftools +munge -C colheaders.tsv -f human_g1k_v37.fasta -s SCZ_2022 -Ob -o PGC3_SCZ.bcf PGC3_SCZ.tsv.gz
```

Liftover VCFs tool:
```
Usage: bcftools +liftover [General Options] -- [Plugin Options]
Options:
   run "bcftools plugin" for a list of common options

Plugin options:
   -s, --src-fasta-ref <file>      source reference sequence in fasta format
   -f, --fasta-ref <file>          destination reference sequence in fasta format
       --set-cache-size <int>      select fasta cache size in bytes
   -c, --chain <file>              UCSC liftOver chain file
       --max-snp-gap <int>         maximum distance to merge contiguous blocks separated by same distance [1]
       --max-indel-gap <int>       maximum distance between contiguous blocks to pad alleles [20]
       --indel-win <int>           maximum distance between two edges of an indel to accept liftover [250]
       --lift-mt                   force liftover of MT/chrMT [automatically determined from contig lengths]
       --no-left-align             do not attempt to left align indels after liftover
       --print-blocks <file>       output contiguous blocks used for the liftOver
       --reject <file>             output variants that cannot be lifted over
   -O, --reject-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
       --write-source              write the source contig/position/alleles for lifted variants

Options for how to update INFO/FORMAT records:
       --flip-tag <string>         INFO annotation flag to record whether alleles are flipped [FLIP]
       --swap-tag <string>         INFO annotation to record when alleles are swapped [SWAP]
       --tags-to-drop <list>       INFO and FORMAT tags to drop when alleles are swapped [INFO/AC,FMT/AC]
       --tags-to-reverse <list>    INFO and FORMAT tags to be reversed when alleles are swapped (must be Number=A,Type=Float)
                                   [INFO/AF:1,FMT/AF:1,FMT/DS:2,FMT/AP1:1,FMT/AP2:1]
       --tags-to-flip <list>       INFO and FORMAT tags that have the sign flipped when alleles are swapped (must be Number=A)
                                   [FMT/EZ,FMT/ES,FMT/ED]
       --tags-genotype <list>      INFO and FORMAT tags with genotype integers like FORMAT/GT (must be Type=Integer)
                                   [INFO/ALLELE_A,INFO/ALLELE_B]

Examples:
      bcftools +liftover -Ob -o output.hg38.bcf input.hg19.bcf -- \
        -s human_g1k_v37.fasta -f Homo_sapiens_assembly38.fasta -c hg19ToHg38.over.chain.gz
      bcftools +liftover -Oz -o chm13v2.0_dbSNPv155.vcf.gz GRCh38_dbSNPv155.vcf.gz -- \
        -s Homo_sapiens_assembly38.fasta -f chm13v2.0.fa -c hg38-chm13v2.over.chain.gz

To obtain UCSC liftOver chain files:
      wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
      wget http://hgdownload.cse.ucsc.edu/goldenPath/hs1/liftOver/hg38-chm13v2.over.chain.gz
```

Best linear unbiased prediction tool:
```
Usage: bcftools +blupx [options] <score.gwas.vcf.gz> [<ldgm.vcf.gz> <ldgm2.vcf.gz> ...]
Plugin options:
   -b, --beta-cov                  frequency-dependent architecture parameter [1e-7]
   -x, --cross-corr                cross ancestry correlation parameter [0.9]
   -a, --alpha-param               alpha parameter [0]
       --tolerance <float>         Tolerance threshold for the conjugate gradient [1e-10]
       --no-jacobi                 Do not use Jacobi preconditioning when solving linear systems with conjugate gradient
       --sample-sizes <list>       List of sample sizes for each input summary statistic [estimated from NS/NC/NE fields]
       --ldgm-vcfs <list>          List of LDGM-VCF files to use
       --ldgm-vcfs-file <file>     File of list of LDGM-VCF files to use
   -e, --exclude EXPR              Exclude sites for which the expression is true (see man page for details)
   -i, --include EXPR              Select sites for which the expression is true (see man page for details)
       --no-version                do not append version and command line to the header
   -o, --output <file>             write output to a file [no output]
   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
   -l, --log <file>                write log to file [standard error]
   -r, --regions <region>          restrict to comma-separated list of regions
   -R, --regions-file <file>       restrict to regions listed in a file
       --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
   -s, --samples <list>            List of summary statitics to include
   -S, --samples-file <file>       File of list of summary statistics to include
   -t, --targets [^]<region>       restrict to comma-separated list of regions. Exclude regions with "^" prefix
   -T, --targets-file [^]<file>    restrict to regions listed in a file. Exclude regions with "^" prefix
       --targets-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]
       --threads <int>             use multithreading with INT worker threads [0]

Examples:
      bcftools +blupx -Ob -o ukb.blup.gwas.bcf -b 2e-7 ukb.gwas.bcf 1kg_ldgm.EUR.bcf
```

Meta-analysis tool:
```
Usage: bcftools +metal [options] <score1.gwas.vcf.gz> <score2.gwas.vcf.gz> [<score3.gwas.vcf.gz> ...]
Plugin options:
       --summaries <file>          list of summary statistics VCFs from file
   -e, --exclude EXPR              Exclude sites for which the expression is true (see man page for details)
   -i, --include EXPR              Select sites for which the expression is true (see man page for details)
       --szw                       perform meta-analysis based on sample-size weighted scheme
                                   rather than inverse-variance weighted scheme
       --het                       perform heterogenity analysis
       --esd                       output effect size direction across studies
       --no-version                do not append version and command line to the header
   -o, --output <file>             write output to a file [no output]
   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
   -r, --regions <region>          restrict to comma-separated list of regions
   -R, --regions-file <file>       restrict to regions listed in a file
       --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
   -t, --targets [^]<region>       restrict to comma-separated list of regions. Exclude regions with "^" prefix
   -T, --targets-file [^]<file>    restrict to regions listed in a file. Exclude regions with "^" prefix
       --targets-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]
       --threads <int>             use multithreading with INT worker threads [0]

Examples:
      bcftools +metal -Ob -o ukb_mvp.gwas.bcf -i ukb.gwas.bcf mvp.gwas.bcf
      bcftools +metal -Ob -o ukb_mvp.gwas.bcf -i 'NS>1000 & AF>0.01 & AF<0.99' ukb.gwas.bcf mvp.gwas.bcf
      bcftools +metal -Ob -o ukb_mvp.gwas.bcf -i 'ID="rs1234" || ID="rs123456" || ID="rs123"' ukb.gwas.bcf mvp.gwas.bcf
```

Installation
============

Install basic tools (Debian/Ubuntu specific if you have admin privileges)
```
sudo apt install wget libcurl4 bcftools r-cran-optparse r-cran-ggplot2 r-cran-data.table
```

Download latest version of [HTSlib](https://github.com/samtools/htslib) and [BCFtools](https://github.com/samtools/bcftools) (if not downloaded already)
```
wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
tar xjvf bcftools-1.16.tar.bz2
```

Download and compile plugins code (make sure you are using gcc version 5 or newer)
```
cd bcftools-1.16/
/bin/rm -f plugins/{score.{c,h},{munge,liftover,metal,blupx}.c}
wget -P plugins https://raw.githubusercontent.com/freeseek/score/master/{score,munge,liftover,blupx,metal}.c
make
/bin/cp bcftools plugins/{munge,liftover,score,metal,blupx}.so $HOME/bin/
wget -P $HOME/bin https://personal.broadinstitute.org/giulio/score/assoc_plot.R
chmod a+x $HOME/bin/assoc_plot.R
```

Make sure the directory with the plugins is available to BCFtools
```
export PATH="$HOME/bin:$PATH"
export BCFTOOLS_PLUGINS="$HOME/bin"
```

Install the GRCh37 human genome reference, cytoband and chain file
```
wget -O- ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz | \
  gzip -d > $HOME/GRCh37/human_g1k_v37.fasta
samtools faidx $HOME/GRCh37/human_g1k_v37.fasta
bwa index $HOME/GRCh37/human_g1k_v37.fasta
wget -P $HOME/GRCh37 http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
wget -P $HOME/GRCh37 http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg19.over.chain.gz
ref="$HOME/GRCh37/human_g1k_v37.fasta"
```

Install the GRCh38 human genome reference (following the suggestion from [Heng Li](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)), cytoband and chain files
```
wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | \
  gzip -d > $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
samtools faidx $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
bwa index $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
wget -P $HOME/GRCh38 http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
wget -P $HOME/GRCh38 http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain.gz
wget -P $HOME/GRCh38 http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
ref="$HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
```

Column Headers Mappings
=======================

Generate column headers mappings from the [MungeSumstats](https://doi.org/10.1093/bioinformatics/btab665) Bioconductor package for importing summary statistics
```
wget https://raw.githubusercontent.com/neurogenomics/MungeSumstats/master/data/sumstatsColHeaders.rda
(Rscript -e 'load("sumstatsColHeaders.rda"); write.table(sumstatsColHeaders, "", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)' | \
  awk -F"\t" -v OFS="\t" '
  ($1~"^ALT" || $1~"^EFF" || $1~"^MINOR" || $1~"^INC" || $1~"T[eE][sS][tT][eE][dD]" || $1=="EA") && $2=="A2" {$2="A1"}
  ($1~"^REF" || $1~"^NON" || $1~"^OTHER" || $1~"^MAJOR" || $1~"^DEC" || $1=="NEA") && $2=="A1" {$2="A2"}
  ($1=="A2FREQ" || $1=="A2FRQ") && $2=="FRQ" {$2="A2FRQ"}
  ($1=="EFFECTIVE_N" || $1=="NEFF") && $2=="N" {$2="NEFF"} {print}'
echo -e "CHR_NAME\tCHR"
echo -e "BP_GRCH38\tBP"
echo -e "CHR_POSITION\tBP"
echo -e "GENPOS\tBP"
echo -e "NAME\tSNP"
echo -e "VARIANT_ID\tSNP"
echo -e "AL1\tA1"
echo -e "AL2\tA2"
echo -e "IMPINFO\tINFO"
echo -e "IMPUTATION\tINFO"
echo -e "R2HAT\tINFO"
echo -e "RSQ\tINFO"
echo -e "EFFECT_WEIGHT\tBETA"
echo -e "INV_VAR_META_BETA\tBETA"
echo -e "ALL_INV_VAR_META_BETA\tBETA"
echo -e "ALL_META_SAMPLE_N\tN"
echo -e "INV_VAR_META_SEBETA\tSE"
echo -e "ALL_INV_VAR_META_SEBETA\tSE"
echo -e "LOG10_P\tLP"
echo -e "LOG10P\tLP"
echo -e "MLOG10P\tLP"
echo -e "P.SE\tP"
echo -e "INV_VAR_META_P\tP"
echo -e "ALL_INV_VAR_META_P\tP"
echo -e "FREQ_EFFECT\tFRQ"
echo -e "ALL_META_AF\tFRQ"
echo -e "NCAS\tN_CAS"
echo -e "NCON\tN_CON"
echo -e "Weight\tNEFF"
echo -e "NEFFDIV2\tNEFFDIV2"
echo -e "HetISq\tHET_I2"
echo -e "HetISqt\tHET_I2"
echo -e "HetPVa\tHET_P"
echo -e "HetPVal\tHET_P"
echo -e "logHetP\tHET_LP"
echo -e "Direction\tDIRE"
echo -e "DIRE\tDIRE") > colheaders.tsv
/bin/rm sumstatsColHeaders.rda
```
Notice that MungeSumstats assigns `A2` rather than `A1` as the effect allele, prompting a correction to revert the mapping to what the original [munge_sumstats.py](https://github.com/bulik/ldsc/blob/master/munge_sumstats.py) had

If your summary statistics file contains headers that cannot be parsed, consider [reporting the issue](https://github.com/neurogenomics/MungeSumstats#future-enhancements) to the MungeSumstats authors

LDGM-VCF Specification
======================

Similar to the [GWAS-VCF specification](https://github.com/MRCIEU/gwas-vcf-specification), an LDGM-VCF file is a VCF file whose header must include the following mandatory INFO fields
```
##INFO=<ID=AA,Number=1,Type=Integer,Description="Ancestral Allele">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=LD_block,Number=1,Type=Integer,Description="Number of LDGM precision matrix">
##INFO=<ID=LD_node,Number=1,Type=Integer,Description="Node corresponding to variant in the LDGM precision matrix">
##INFO=<ID=LD_diagonal,Number=1,Type=Float,Description="Weight of the node in the LDGM precision matrix">
##INFO=<ID=LD_neighbors,Number=.,Type=Integer,Description="Nodes of the neighbors in the LDGM precision matrix">
##INFO=<ID=LD_weights,Number=.,Type=Float,Description="Weights of the edges in the LDGM precision matrix">
```

There should be only one alternate allele per line and the `AA` field must be a number equal to 0 if the ancestral allele is the reference allele and 1 if the ancestral allele is the alternate allele. The `LD_block` field must be a non-negative integer monotonically increasing across variants and indicating which LDGM matrix a given variant is part of. The `LD_node` field must be a non-negative integer indicating which node of the LDGM matrix a variant corresponds to. It is allowed for variants in perfect linkage disequilibrium to have the same `LD_block` and `LD_node` values. The `LD_node` numbers across variants do not need to be monotonically increasing and it is okay for some `LD_node` numbers to be missing from a given LDGM matrix. The `LD_diagonal` must be a number equal or larger then one. The `LD_neighbors` and `LD_weigths` arrays must have the same length. The integer numbers within the `LD_neighbors` arrays must all greater than the `LD_node` number, as the LDGM matrix, given its symmetry, must be stored in triangular upper format to save space. The floating point numbers within the `LD_weigths` arrays must be non-zero

| #CHROM | POS   | ID          | REF | ALT | QUAL | FILTER | INFO                                                                                                                               |
|--------|-------|-------------|-----|-----|------|--------|------------------------------------------------------------------------------------------------------------------------------------|
| chr1   | 16719 | rs62636367  | T   | A   | .    | .      | AA=0;AF=0.0626;LD_block=0;LD_node=4;LD_diagonal=1.55379;LD_neighbors=6,12,21,52;LD_weights=-0.319217,-0.466229,-0.066764,-0.247807 |
| chr1   | 16841 | rs62636368  | G   | T   | .    | .      | AA=0;AF=0.0855;LD_block=0;LD_node=6;LD_diagonal=1.73014;LD_neighbors=12;LD_weights=-0.914626                                       |
| chr1   | 16856 | rs3891260   | A   | G   | .    | .      | AA=0;AF=0.0308;LD_block=0;LD_node=7;LD_diagonal=1                                                                                  |
| chr1   | 16949 | rs199745162 | A   | C   | .    | .      | AA=0;AF=0.3668;LD_block=0;LD_node=8;LD_diagonal=3.26079;LD_neighbors=10,18,57,114;LD_weights=-1.6973,-1.10987,-0.135282,-0.048439  |
| chr1   | 17005 | rs201833382 | A   | G   | .    | .      | AA=0;AF=0.0656;LD_block=0;LD_node=9;LD_diagonal=1.14963;LD_neighbors=35,94,5358;LD_weights=-0.332079,-0.185058,-0.1273             |

Representing the ancestral allele with a number rather than with a string referring to the ancestral allele as done by the [International Genome Sample Resource](https://www.internationalgenome.org/wiki/Analysis/vcf4.0/) is helpful both to improve processing speed and for compatibility with the operation of left-aligning indels that can be performed with the command `bcftools norm --fasta-ref`

Variants in perfect linkage disequilibrium with the same `LD_block` and `LD_node` values must also have the same `LD_neighbors` and `LD_weights` array values, while they can have different `AA` values. This will cause a slight loss of redundancy as approximately 15% of variants can be considered redundant due to perfect linkage disequilibrium. The signs of the weights of the LDGM matrix refer to the derived alleles, which in approximately 85% of cases is the alternate allele

The `ID` field does not need to be filled as matrices from and LDGM-VCF file and summary statistics from a GWAS-VCF file will be unequivocally matched using genomic position, reference and alternate alleles

LDGM Matrices
=============

Linkage disequilibrium graphical models (LDGM) precision matrices for [1,361 intervals](https://github.com/jmacdon/LDblocks_GRCh38) computed for the GRCh38 genome can be downloaded from [here](https://ldgm.readthedocs.io/en/latest/introduction.html). However, SNP list files are provided without position information, so we need to first recover this information to be able to match the SNPs to the SNPs in a summary statistics file following the GWAS-VCF specification. You can download the LDGM-VCF precision matrices [here](http://software.broadinstitute.org/software/score)

The following code will generate updated SNP list files with recovered position information and knowledge of whether the ancestral allele was the reference or the alternate allele by tracing back the steps used to generate the provided SNP lists from the [LDGM paper](https://github.com/awohns/ldgm_paper)
```
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/BED/bed_chr_{1..22}.bed.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{1..22}.filtered.shapeit2-duohmm-phased.vcf.gz{,.tbi}
wget -O snplist.tar.gz https://www.dropbox.com/sh/raw/1huaxgad2bjjv9a/AAD9YEljtU3TxYum3qPxJIp6a/ldgm/snplist.tar.gz?dl=0
tar xzvf snplist.tar.gz

mkdir -p ids
for chr in {1..22}; do zcat bed_chr_$chr.bed.gz | tail -n+2 | cut -f3,4 | sort -k2,2 > bed_chr_$chr.tsv; done
for chr in {1..22}; do
  for file in snplist/1kg_chr${chr}_[0-9]*_[0-9]*.snplist; do
    lbl=${file%.snplist};
    lbl=${lbl#*1kg_};
    cut -d, -f9 $file | sort | join -1 1 -2 2 - bed_chr_$chr.tsv | tr ' ' ',' > ids/$lbl.csv
  done
done
/bin/rm bed_chr_{1..22}.tsv

mkdir -p afs
inc="AC_EUR_unrel/AN_EUR_unrel>.01 && AC_EUR_unrel/AN_EUR_unrel<=.99 || AC_EAS_unrel/AN_EAS_unrel>=.01 && AC_EAS_unrel/AN_EAS_unrel<=.99 || AC_AMR_unrel/AN_AMR_unrel>=.01 && AC_AMR_unrel/AN_AMR_unrel<=.99 || AC_SAS_unrel/AN_SAS_unrel>=.01 && AC_SAS_unrel/AN_SAS_unrel<=.99 || AC_AFR_unrel/AN_AFR_unrel>=.01 && AC_AFR_unrel/AN_AFR_unrel<=.99"
fmt="%REF,%ALT,%AC_EUR_unrel,%AN_EUR_unrel,%AC_EAS_unrel,%AN_EAS_unrel,%AC_AMR_unrel,%AN_AMR_unrel,%AC_SAS_unrel,%AN_SAS_unrel,%AC_AFR_unrel,%AN_AFR_unrel,%POS\n"
for chr in {1..22}; do
  vcf="CCDG_14151_B01_GRM_WGS_2020-08-05_chr$chr.filtered.shapeit2-duohmm-phased.vcf.gz"
  for file in snplist/1kg_chr${chr}_[0-9]*_[0-9]*.snplist; do
    lbl=${file%.snplist};
    lbl=${lbl#*1kg_};
    reg=${lbl/_/:};
    reg=${reg/_/-};
    bcftools query -f "$fmt" -i "$inc" -r $reg $vcf | \
      awk -F, '{printf "%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,+,%d\n",$1,$2,$3/$4,$5/$6,$7/$8,$9/$10,$11/$12,$13;
        printf "%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,-,%d\n",$2,$1,($4-$3)/$4,($6-$5)/$6,($8-$7)/$8,($10-$9)/$10,($12-$11)/$12,$13}' | \
      sed 's/-0/0/g;s/0,/,/g;s/0,/,/g;s/0,/,/g' > afs/$lbl.csv
  done
done

mkdir -p out
for file in snplist/1kg_chr[0-9]*_[0-9]*_[0-9]*.snplist; do
  lbl=${file%.snplist};
  lbl=${lbl#*1kg_};
  awk -F, 'BEGIN {x["site_ids"]="position"; x["NA"]="NA"}
    NR==FNR {x[$1]=$2} NR>FNR {print $0","x[$9]}' ids/$lbl.csv $file | \
  awk -F, -v OFS=, 'BEGIN {y["anc_alleles,deriv_alleles,EUR,EAS,AMR,SAS,AFR,position"]="swap"; last=0}
    NR==FNR {str=$1","$2","$3","$4","$5","$6","$7; if (str in x) x[str]=x[str]","$9; else x[str]=$9; y[str","$9]=$8}
    NR>FNR {str=$2","$3","$4","$5","$6","$7","$8; if ($10=="NA" && str in x) {
    split(x[str],a,","); for (i=1; i<=length(a); i++) if (a[i]>last) {$10=a[i]; break}}
    $11=y[str","$10]; print; last=$10}' afs/$lbl.csv - > out/1kg_$lbl.snplist
done

/bin/rm -r snplist ids afs
```

With the updated SNP list files we can format the LDGM precision matrices into LDGM-VCF files
```
wget -O AFR.tar.gz https://www.dropbox.com/sh/raw/1huaxgad2bjjv9a/AADu-h_GZF7FI2FoNJYN9t9Oa/ldgm/AFR.tar.gz?dl=0
wget -O AMR.tar.gz https://www.dropbox.com/sh/raw/1huaxgad2bjjv9a/AADhcJm-THCOX5gpCKqZmvpva/ldgm/AMR.tar.gz?dl=0
wget -O EAS.tar.gz https://www.dropbox.com/sh/raw/1huaxgad2bjjv9a/AADCBA9TrjQoSJiF4fbJ2oLZa/ldgm/EAS.tar.gz?dl=0
wget -O EUR.tar.gz https://www.dropbox.com/sh/raw/1huaxgad2bjjv9a/AAB8i85pOY-XVNPnQ9NUwUaAa/ldgm/EUR.tar.gz?dl=0
wget -O SAS.tar.gz https://www.dropbox.com/sh/raw/1huaxgad2bjjv9a/AADbbgk0VErJ_dXC7D1L-p3ga/ldgm/SAS.tar.gz?dl=0

(echo "##fileformat=VCFv4.2"
echo "##INFO=<ID=AA,Number=1,Type=Integer,Description=\"Ancestral Allele\">"
echo "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"
echo "##INFO=<ID=LD_block,Number=1,Type=Integer,Description=\"Number of LDGM precision matrix\">"
echo "##INFO=<ID=LD_node,Number=1,Type=Integer,Description=\"Node corresponding to variant in the LDGM precision matrix\">"
echo "##INFO=<ID=LD_diagonal,Number=1,Type=Float,Description=\"Weight of the node in the LDGM precision matrix\">"
echo "##INFO=<ID=LD_neighbors,Number=.,Type=Integer,Description=\"Nodes of the neighbors in the LDGM precision matrix\">"
echo "##INFO=<ID=LD_weights,Number=.,Type=Float,Description=\"Weights of the edges in the LDGM precision matrix\">"
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO") > tmp.vcf
for anc in AFR AMR EAS EUR SAS; do
  tar xzvf $anc.tar.gz
  (bcftools reheader --fai $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai --temp-prefix ./bcftools. tmp.vcf
  ls out/1kg_chr[0-9]*_[0-9]*_[0-9]*.snplist | \
    sed 's/out\/1kg_chr//' | \
    sort -t_ -k1,1n -k2,2n | \
    sed 's/^\([0-9]*\)\(_[0-9]*_[0-9]*\)\.snplist$/chr\1 out\/1kg_chr\1\2.snplist '$anc'\/1kg_chr\1\2.'$anc'.edgelist/' | \
    cat -n | \
  while read block chr snpfile edgefile; do
    awk -F, -v anc=$anc -v chr=$chr -v block=$((block-1)) '
      NR==FNR && $1==$2 {x[$1]=$3} NR==FNR && $1!=$2 {y[$1]=y[$1]" "$2; z[$1]=z[$1]" "$3}
      NR>FNR && FNR==1 {for (i=1; i<=NF; i++) f[$i] = i}
      NR>FNR && FNR>1 && ($1 in x || $1 in y) {ref=$(f["anc_alleles"]); alt=$(f["deriv_alleles"]);
      pos=$(f["position"]); aa=0; af=$(f[anc]); node=$(f["index"]); score=x[node];
      if ($(f["swap"])=="-") {ref=$(f["deriv_alleles"]); alt=$(f["anc_alleles"]); aa=1; af=1-af}
      printf "%s\t%d\t.\t%s\t%s\t.\t.\tAA=%d;AF=%f;LD_block=%d;LD_node=%d;LD_diagonal=%s",chr,pos,ref,alt,aa,af,block,node,score
      if ($1 in y) {neighbors=substr(y[$1],2); gsub(" ", ",", neighbors);
      weights=substr(z[$1],2); gsub(" ", ",", weights); printf ";LD_neighbors=%s;LD_weights=%s",neighbors,weights}
      printf "\n"}' $edgefile $snpfile
  done) | bcftools view --no-version -Ob | \
  tee 1kg_ldgm.$anc.bcf | \
  bcftools index --force --output 1kg_ldgm.$anc.bcf.csi
  /bin/rm -r $anc
done
/bin/rm tmp.vcf
```

You can recover the LDGM matrix in the original format compatible with the LDGM [readedgelist](https://github.com/awohns/ldgm/blob/main/matlab/utility/readedgelist.m) function
```
bcftools query -i "LD_block=135" -f "%LD_node\t%LD_diagonal\t%LD_neighbors\t%LD_weights\n" -r chr2:55438332-59565357 1kg_ldgm.EUR.bcf | \
  awk -F"\t" -v OFS=, '{print $1,$1,$2} $3!="." {split($3,a,","); split($4,b,","); for (i=1; i<=length(a); i++) print $1,a[i],b[i]}' | \
  sort -t, -k1,1n -k2,2n | uniq
```

To split the LDGM-VCF file in 1,361 LDGM-VCF files containing each block separately
```
wget -O- https://raw.githubusercontent.com/jmacdon/LDblocks_GRCh38/master/data/EUR_LD_blocks.bed | \
  awk 'NR>1 {printf "%s:%d:%d\n",$1,$2,$3}' > EUR_LD_blocks.txt
ulimit -n 2048
bcftools +scatter --no-version -Ob 1kg_ldgm.EUR.bcf -o EUR -S EUR_LD_blocks.txt -o LD_blocks/ -p 1kg_ldgm.EUR.
```

Compute polygenic scores
========================

The BCFtools score plugin can input summary statistics files in a variety of formats, including those following the [GWAS-VCF specification](https://github.com/MRCIEU/gwas-vcf-specification), those following the [GWAS-SSF specification](https://ebispot.github.io/gwas-blog/new-standard-for-gwas-summary-statistics), and more in general most summary statistics files formatted as text tables with a header indicating which column to use. For GWAS-SSF and table summary statistiscs files, BCFtools score will automatically recognize the columns and attempt to match variants by chromosome and position if available and then by marker name if the genomic position is unavailable in the summary statistics file. Multiple summary statistics files can be input at once except you cannot mix GWAS-VCF summary statistics files with other files. If multiple summary statistics are present in a GWAS-VCF, all will be scored independently

One advantage of the BCFtools score plugin is that it can be readily used on imputation VCFs without further format conversion. It will work with Minimac3, Minimac4, Beagle5, and IMPUTE5 output VCFs and more in general with any VCF including any of the following format fields

| FORMAT  | Description                                      |
|---------|--------------------------------------------------|
| AP1/AP2 | ALT allele probability of first/second haplotype |
| HDS     | Estimated Haploid Alternate Allele Dosage        |
| GP      | Estimated Genotype Probability                   |
| DS      | Genotype dosage                                  |
| GT      | Genotype                                         |

Convert summary statistics
==========================

The BCFtools munge plugin, inspired by the [MungeSumstats](https://doi.org/10.1093/bioinformatics/btab665) tool from Alan Murphy which is itself inspired by the `munge_sumstats.py` script in [ldsc](https://github.com/bulik/ldsc) from Brendan Bulik-Sullivan, allows the majority of summary statitsics files available to the scientific community to be converted to summary statistics files following the [GWAS-VCF specification](https://github.com/MRCIEU/gwas-vcf-specification)

While being an alternative to MungeSumStats and `munge_sumstats.py`, the BCFtools munge plugin does not support the same number of features with some differences highlighted in the following table

| Feature                               | MungeSumStats | BCFtools +munge |
|---------------------------------------|---------------|-----------------|
| outputs GWAS-VCF                      | YES           | YES             |
| handles either tab or space delimited | YES           | YES             |
| handles header name synonyms          | YES           | YES             |
| remove strand-ambiguous SNPs          | YES           | NO              |
| check for allele flipping from AF     | YES           | NO              |
| check whether A1 or A2 is reference   | NO            | YES             |
| assumes as effect allele ...          | A2            | A1              |

Notice however that for many indels it is impossible to retrieve which allele is the reference allele if the table does not explicitly specify which allele is the reference allele as sometimes both alleles can match the reference sequence, a problem that the VCF specification was designed to solve

To convert a given summary statistics file generated by [PLINK](https://www.cog-genomics.org/plink/) you can simply run a command like the following
```
wget https://raw.githubusercontent.com/neurogenomics/MungeSumstats/master/inst/extdata/ieu-a-298.tsv.gz
bcftools +munge --no-version -c PLINK -f $HOME/GRCh37/human_g1k_v37.fasta -s ieu-a-298 ieu-a-298.tsv.gz
```

If you want to convert to a different reference genome
```
zcat ieu-a-298.tsv.gz | \
bcftools +munge --no-version -Ou -c PLINK --fai $HOME/GRCh37/human_g1k_v37.fasta.fai -s ieu-a-298 |
bcftools +liftover --no-version -Ob -- \
  -f $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -c $HOME/GRCh38/hg19ToHg38.over.chain.gz | \
tee ieu-a-298.hg38.bcf | \
bcftools index --force --output ieu-a-298.hg38.bcf.csi
```

For summary statistics files following a less specific column header format, you can use a comprehensive column headers mapping
```
wget https://raw.githubusercontent.com/neurogenomics/MungeSumstats/master/inst/extdata/eduAttainOkbay.txt
bcftools +munge --no-version -Ou -C colheaders.tsv --fai $HOME/GRCh37/human_g1k_v37.fasta.fai -s eduAttain eduAttainOkbay.txt | \
bcftools +liftover --no-version -Ou -- \
  -f $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -c $HOME/GRCh38/hg19ToHg38.over.chain.gz | \
bcftools sort -Ob | \
tee eduAttainOkbay.hg38.bcf | \
bcftools index --force --output eduAttainOkbay.hg38.bcf.csi
```

For summary statistics files including indels, you will need to provide both references when performing the liftover
```
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_10k_SNPs.zip
unzip -p COVID19_HGI_10k_SNPs.zip COVID19_HGI_A2_ALL_20210107.10k.b37.txt.gz | \
bcftools +munge --no-version -Ou -C colheaders.tsv --fai $HOME/GRCh37/human_g1k_v37.fasta.fai -s COVID_2021 | \
bcftools +liftover --no-version -Ou -- \
  -s $HOME/GRCh37/human_g1k_v37.fasta \
  -f $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -c $HOME/GRCh38/hg19ToHg38.over.chain.gz | \
bcftools sort -Ob | \
tee COVID19_HGI_A2_ALL_20210107.10k.hg38.bcf | \
bcftools index --force --output COVID19_HGI_A2_ALL_20210107.10k.hg38.bcf.csi
```

Liftover VCFs
=============

The BCFtools liftover plugin is inspired by the Picard [LiftoverVcf](https://gatk.broadinstitute.org/hc/en-us/articles/360037060932-LiftoverVcf-Picard-) tool, written by Alec Wysoker, Benjamin Bimber, Tim Fennell, and Yossi Farjoun, and allows to liftover VCFs from one reference to another including summary statistics files following the [GWAS-VCF specification](https://github.com/MRCIEU/gwas-vcf-specification). Beyond being much faster than the Picard LiftoverVcf tool, the BCFtools liftover plugin supports several additional features summarized in the following table

| Feature                            | Picard LiftoverVcf | BCFtools +liftover |
|------------------------------------|--------------------|--------------------|
| SNPs                               | YES                | YES                |
| indels                             | YES                | YES                |
| left align indels after liftover   | YES                | YES                |
| sort records after liftover        | YES                | NO                 |
| SNPs at 1bp gaps in the chain file | NO                 | YES                |
| flips alleles when changing strand | YES                | YES                |
| swaps SNP alleles when needed      | YES                | YES                |
| swaps indel alleles when needed    | NO                 | YES                |
| adds reference alleles when needed | NO                 | YES                |
| handles GT/PL/AD records           | bi-allelic only    | YES                |
| handles Number=G/Number=R records  | NO                 | YES                |
| reverses Number=A records          | only AF-like       | YES                |
| handles EZ/ES/AF GWAS-VCF records  | NO                 | YES                |
| flexible with contig names         | NO                 | YES                |
| can input a VCF as a file stream   | NO                 | YES                |
| can input and output binary VCFs   | NO                 | YES                |
| loads whole reference in memory    | YES                | NO                 |

At the time of this writing the BCFtools liftover plugin is the only liftover tool that handles indels at short tandem repeats correctly even in the case the two reference genomes represent different alleles as well as SNPs that fall within 1bp gaps between contiguous blocks from the same chain. When applied it to variants from the 1000 Genomes project phase 3 (GRCh37 sites file available [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)), the BCFtools liftover plugin and Picard LiftoverVcf (run with option `--LIFTOVER_MIN_MATCH 0.0` and `--RECOVER_SWAPPED_REF_ALT true`) have the following statistics

|                    | Feature   | SNP        | INDEL     | SNP,INDEL | MNP |
|--------------------|-----------|------------|-----------|-----------|-----|
| Tool               | Total     | 81,377,202 | 3,299,133 | 65,871    | 123 |
|--------------------|-----------|------------|-----------|-----------|-----|
|                    | Rejected  | 24,130     | 1,113     | 34        | 0   |
| BCFtools +liftover | Ref added | 162        | 2,289     | 96        | 0   |
|                    | Swapped   | 36,967     | 5,904     | 123       | 2   |
|--------------------|-----------|------------|-----------|-----------|-----|
|                    | Rejected  | 29,472     | 4,890     | 144       | 2   |
| Picard LiftoverVcf | Incorect  | 0          | 4,360     | 118       | 10  |
|                    | Swapped   | 31,787     | 0         | 0         | 0   |

To be able to swap reference and alternate alleles for indels when needed, the BCFtools liftover plugin uses the source reference to first extend all the alleles until they have a unique representation that makes it mathematically impossible to match the wrong allele after liftover to the destination reference. To further recover more indels, if one edge of the sequence being lifted over falls within any of the contiguous blocks from one of the chains but not the other end, the BCFtools liftover plugin will further extend the sequence until both ends fall within contiguous blocks from the same chain

It can be tested as follows
```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz{,.tbi}
bcftools +liftover --no-version -Ou ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.bcf -- \
  -s $HOME/GRCh37/human_g1k_v37.fasta \
  -f $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -c $HOME/GRCh38/hg19ToHg38.over.chain.gz \
  --reject ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.reject.bcf \
  --reject-type b \
  --write-src | \
bcftools sort -Ob | \
tee ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.hg38.bcf | \
bcftools index --force --output ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.hg38.bcf.csi
```

Compute best linear unbiased predictor
======================================

The BCFtools blupx plugin is inspired by the [BLUPx-ldgm](https://doi.org/10.1101/2022.09.06.506858) software, written and designed by Pouria Salehi Nowbandegani, Anthony Wilder Wohns, and Luke O’Connor, and it will apply the [best linear unbiased prediction](https://en.wikipedia.org/wiki/Best_linear_unbiased_prediction) (BLUP) model to compute improved polygenic weights starting from summary statistics following the [GWAS-VCF specification](https://github.com/MRCIEU/gwas-vcf-specification) following the MATLAB code from the [LDGM](https://github.com/awohns/ldgm) repository

First of all, compute the number of markers shared between your LDGM-VCF file and your GWAS-VCF file:
```
b=1e-7
bcftools isec -n 2 <score.gwas.vcf.gz> <ldgm.vcf.gz> -w 2 | bcftools query -f "%AF\n" | awk -v b=$b '{S+=2*$1*(1-$1)} END {print S*b}'
```
The output will give you the expected heritability from common SNPs. You can adjust the value `b` to yield the expected heritabilty for the trait you are working with. Once you have the correct value `b` you can then generate the BLUP loadings
```
bcftools +blupx \
  --no-version \
  --beta-cov $b \
  <score.gwas.vcf.gz> \
  <ldgm.vcf.gz> \
  --output-type b \
  --output <score.blup$b.vcf.gz>
```
You can also generate BLUP loadings for different values of `b` and then merge the output GWAS-VCFs into a single GWAS-VCF file that you can then use to compare the performance of different choices for `b`

Run meta-analysis
=================

The BCFtools metal plugin is inspired by the [METAL](http://doi.org/10.1093/bioinformatics/btq340) software written by Goncalo Abecasis and it performs fixed effect meta-analyses from summary statistics following the [GWAS-VCF specification](https://github.com/MRCIEU/gwas-vcf-specification) using either the inverse-variance weighted (IVW) scheme or the sample-size weighted (SZW) scheme. Both softwares can filter variants, METAL through [filtering conditions](https://genome.sph.umich.edu/wiki/METAL_Documentation#Filtering) and BCFtools metal through [filtering expressions](http://samtools.github.io/bcftools/howtos/filtering.html). There are a few differences between the two approaches though, summarized in the following table

| Feature                          | METAL   | BCFtools +metal |
|----------------------------------|---------|-----------------|
| inverse-variance weighted scheme | YES     | YES             |
| sample-size weighted scheme      | YES     | YES             |
| heterozygosity test              | YES     | YES             |
| filter variants                  | YES     | YES             |
| genomic control                  | YES     | NO              |
| corrects for samples overlap     | YES     | NO              |
| match variants by ID             | YES     | NO              |
| match variants by position       | NO      | YES             |
| computes N_eff for binary traits | NO      | YES             |
| input and output GWAS-VCF        | NO      | YES             |
| input p-values in log space      | NO      | YES             |
| output p-values in log space     | YES     | YES             |
| output FreqSE/MinFreq/MaxFreq    | YES     | NO              |
| output HetChiSq/HetDf            | YES     | NO              |
| output NS/NC/AC                  | NO      | YES             |
| output sorted variants           | NO      | YES             |
| multiple phenotypes at once      | NO      | YES             |

If some missing features are important to you, contact the [author](mailto:giulio.genovese@gmail.com) to discuss adding option to the BCFtools metal plugin. The latter is meant to function as a simplified version of the original METAL software allowing to perform the most common meta-analyses while inputting and outputting files in a standardized file format. It requires summary statistics to be properly formatted which can be accomplished using `bcftools +munge` and `bcftools +liftover`

This is an example to compare how to use the original METAL software and the BCFtools metal plugin
```
wget http://csg.sph.umich.edu/abecasis/metal/download/GlucoseExample-original.tar.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg17/liftOver/hg17ToHg38.over.chain.gz
tar xzvf GlucoseExample-original.tar.gz
cd GlucoseExample/
echo -e "chr2\t243018229\t0\t0\t0\nchr7\t158628139\t0\t0\t0\nchr11\t134452384\t0\t0\t0" > hg17.fai
bcftools +munge --no-version -Ou -C colheaders.tsv --fai hg17.fai -s glucose DGI_three_regions.txt | \
  bcftools +liftover --no-version -Ou -- \
    -f $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg17ToHg38.over.chain.gz | \
  bcftools sort -Ob -o DGI_three_regions.bcf
zcat MAGIC_FUSION_Results.txt.gz | sed '1s/FREQ_EFFECT/FREQ_EFFECT_ALLELE/;s/GEN/1.0/' | \
  bcftools +munge --no-version -Ou -C colheaders.tsv --fai hg17.fai -s glucose | \
  bcftools +liftover --no-version -Ou -- \
    -f $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg17ToHg38.over.chain.gz | \
  bcftools sort -Ob -o MAGIC_FUSION_Results.bcf
cat magic_SARDINIA.tbl | sed '1s/AL1/A1/;1s/AL2/A2/' | \
  bcftools +munge --no-version -Ou -C colheaders.tsv --fai hg17.fai -s glucose --ns 4108 | \
  bcftools +liftover --no-version -Ou -- \
    -f $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg17ToHg38.over.chain.gz | \
  bcftools sort -Ob -o magic_SARDINIA.bcf
echo DGI_three_regions.bcf MAGIC_FUSION_Results.bcf magic_SARDINIA.bcf | xargs -n1 bcftools index --force
```

Run the inverse-variance weighted meta-analysis
```
$ sed -i 's/\r$//;s/^SCHEME   SAMPLESIZE$/# SCHEME   STDERR/;s/^# AVERAGEFREQ ON$/AVERAGEFREQ ON/;s/^ANALYZE$/LOGPVALUE ON\nANALYZE HETEROGENEITY/' metal.txt
$ metal < metal.txt
$ awk 'NR==1 || $8<-7.3' METAANALYSIS1.TBL | column -t
MarkerName  Allele1  Allele2  Freq1   FreqSE  Effect   StdErr  log(P)  Direction  HetISq  HetChiSq  HetDf  logHetP
rs853781    a        g        0.5160  0.0385  -0.1061  0.0192  -7.46   ---        0.0     0.801     2      -0.1739
rs853789    a        g        0.3810  0.0317  -0.1245  0.0200  -9.30   ---        65.3    5.771     2      -1.253
rs537183    t        c        0.5982  0.0536  0.1128   0.0201  -7.73   +++        47.2    3.785     2      -0.8219
rs853787    t        g        0.6184  0.0314  0.1244   0.0201  -9.21   +++        65.8    5.849     2      -1.27
rs560887    t        c        0.3406  0.0345  -0.1359  0.0203  -10.69  ---        74.1    7.712     2      -1.675
rs475612    t        c        0.4014  0.0498  -0.1107  0.0200  -7.49   ---        47.0    3.775     2      -0.8197
rs853773    a        g        0.5037  0.0187  -0.1152  0.0199  -8.16   ---        54.0    4.347     2      -0.944
rs502570    a        g        0.4020  0.0534  -0.1131  0.0201  -7.76   ---        47.5    3.809     2      -0.8272
rs557462    t        c        0.5977  0.0531  0.1126   0.0201  -7.70   +++        46.9    3.767     2      -0.818
rs563694    a        c        0.5977  0.0541  0.1122   0.0201  -7.66   +++        46.4    3.732     2      -0.8103
$ bcftools +metal --het --esd DGI_three_regions.bcf MAGIC_FUSION_Results.bcf magic_SARDINIA.bcf | \
  bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%NS\t%ES\t%SE\t%LP\t%AF\t%I2\t%CQ\t%ED]\n" -i 'LP>7.3' -H | \
  sed 's/^# //;s/\[[0-9]*\]//g' | column -t
CHROM  POS        ID        REF  ALT  glucose:NS  glucose:ES  glucose:SE  glucose:LP  glucose:AF  glucose:I2  glucose:CQ  glucose:ED
chr2   168906638  rs560887  T    C    6796        0.135859    0.0202669   10.6914     0.659419    74.0672     1.67469     +++
chr2   168917561  rs563694  C    A    6796        0.112246    0.0200544   7.6616      0.59773     46.4042     0.810315    +++
chr2   168918136  rs537183  C    T    6796        0.112759    0.0200544   7.72579     0.598243    47.1627     0.821946    +++
chr2   168918449  rs502570  A    G    6796        0.113055    0.0200544   7.76289     0.598025    47.4977     0.827192    +++
chr2   168920236  rs475612  T    C    6796        0.110729    0.0200347   7.48671     0.59857     47.0206     0.819742    +++
chr2   168921085  rs557462  C    T    6796        0.112572    0.0200558   7.70135     0.597707    46.9099     0.818032    +++
chr2   168944978  rs853789  A    G    6796        0.124542    0.0200248   9.30184     0.619019    65.343      1.25312     +++
chr2   168945742  rs853787  G    T    6796        0.124445    0.0201217   9.20578     0.618385    65.8059     1.27009     +++
chr2   168949811  rs853781  A    G    6796        0.106109    0.0192413   7.45665     0.48399     0           0.173932    +++
chr2   168957837  rs853773  A    G    6796        0.115158    0.0198751   8.16308     0.496252    53.995      0.944016    +++
```

Run the sample-size weighted meta-analysis
```
$ sed -i 's/\r$//;s/^# SCHEME   STDERR$/SCHEME   SAMPLESIZE/;s/^# AVERAGEFREQ ON$/AVERAGEFREQ ON/;s/^ANALYZE$/LOGPVALUE ON\nANALYZE HETEROGENEITY/' metal.txt
$ metal < metal.txt
$ awk 'NR==1 || $8<-7.3' METAANALYSIS1.TBL | column -t
MarkerName  Allele1  Allele2  Freq1   FreqSE  Weight   Zscore  log(P)  Direction  HetISq  HetChiSq  HetDf  logHetP
rs853781    a        g        0.5229  0.0375  6796.00  -5.532  -7.50   ---        0.0     0.156     2      -0.0339
rs853789    a        g        0.3869  0.0310  6796.00  -6.395  -9.79   ---        49.3    3.946     2      -0.8569
rs537183    t        c        0.5884  0.0524  6796.00  5.726   -7.99   +++        31.6    2.923     2      -0.6348
rs853787    t        g        0.6128  0.0307  6796.00  6.401   -9.81   +++        51.0    4.082     2      -0.8864
rs569805    a        t        0.4207  0.0573  6796.00  -5.509  -7.44   ---        24.0    2.630     2      -0.5712
rs560887    t        c        0.3462  0.0336  6796.00  -6.853  -11.14  ---        62.9    5.392     2      -1.171
rs475612    t        c        0.4107  0.0487  6796.00  -5.604  -7.68   ---        29.5    2.835     2      -0.6157
rs579060    t        g        0.5793  0.0573  6796.00  5.506   -7.44   +++        23.9    2.629     2      -0.5708
rs853773    a        g        0.5066  0.0176  6796.00  -5.849  -8.31   ---        19.7    2.490     2      -0.5407
rs508506    a        c        0.4207  0.0573  6796.00  -5.464  -7.33   ---        22.9    2.593     2      -0.563
rs502570    a        g        0.4118  0.0522  6796.00  -5.720  -7.97   ---        30.8    2.889     2      -0.6273
rs552976    a        g        0.4222  0.0579  6796.00  -5.543  -7.53   ---        29.1    2.822     2      -0.6127
rs557462    t        c        0.5880  0.0519  6796.00  5.724   -7.98   +++        29.9    2.854     2      -0.6198
rs486981    a        g        0.4207  0.0573  6796.00  -5.493  -7.40   ---        23.4    2.612     2      -0.5671
rs563694    a        c        0.5878  0.0529  6796.00  5.694   -7.91   +++        31.5    2.919     2      -0.6338
$ bcftools +metal --szw --het --esd DGI_three_regions.bcf MAGIC_FUSION_Results.bcf magic_SARDINIA.bcf | \
  bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%EZ\t%LP\t%AF\t%NE\t%I2\t%CQ\t%ED]\n" -i 'LP>7.3' -H | \
  sed 's/^# //;s/\[[0-9]*\]//g' | column -t
CHROM  POS        ID        REF  ALT  glucose:EZ  glucose:LP  glucose:AF  glucose:NE  glucose:I2  glucose:CQ  glucose:ED
chr2   168906638  rs560887  T    C    6.85292     11.1405     0.653848    6796        62.905      1.17076     +++
chr2   168917561  rs563694  C    A    5.69396     7.90614     0.587813    6796        31.4791     0.633814    +++
chr2   168918136  rs537183  C    T    5.72613     7.98824     0.588417    6796        31.5878     0.63482     +++
chr2   168918449  rs502570  A    G    5.72025     7.97318     0.588236    6796        30.7689     0.627312    +++
chr2   168920236  rs475612  T    C    5.60442     7.67995     0.589339    6796        29.4595     0.615667    +++
chr2   168921085  rs557462  C    T    5.72385     7.98239     0.587981    6796        29.9276     0.619779    +++
chr2   168925639  rs486981  A    G    5.49323     7.4038      0.579313    6796        23.4172     0.567091    +++
chr2   168926370  rs569805  A    T    5.50932     7.44343     0.579313    6796        23.9638     0.571168    +++
chr2   168926529  rs579060  G    T    5.50649     7.43645     0.579313    6796        23.9183     0.570826    +++
chr2   168928445  rs508506  A    C    5.46436     7.33295     0.579313    6796        22.8653     0.563034    +++
chr2   168934928  rs552976  A    G    5.54344     7.52786     0.57781     6796        29.1193     0.612712    +++
chr2   168944978  rs853789  A    G    6.39475     9.79368     0.613109    6796        49.3162     0.856871    +++
chr2   168945742  rs853787  G    T    6.4013      9.81231     0.612819    6796        51.0051     0.886407    +++
chr2   168949811  rs853781  A    G    5.53209     7.4997      0.477075    6796        0           0.0338982   +++
chr2   168957837  rs853773  A    G    5.84872     8.30508     0.493419    6796        19.6833     0.540728    +++
```

Plot results for each study individually and for the inverse-variance weighted meta-analysis meta-analysis
```
bcftools +metal --no-version --het DGI_three_regions.bcf MAGIC_FUSION_Results.bcf magic_SARDINIA.bcf -Ob | \
  tee METAANALYSIS1.bcf | bcftools index --force --output METAANALYSIS1.bcf.csi
for pfx in magic_SARDINIA DGI_three_regions MAGIC_FUSION_Results METAANALYSIS1; do
  for reg in chr2:168411820-169393292 chr7:43699132-44694724 chr11:92476687-93473731; do
    assoc_plot.R --cytoband $HOME/GRCh38/cytoBand.txt.gz --vcf $pfx.bcf --region $reg --png $pfx.${reg%:[0-9]*-[0-9]*}.png
  done
done
```

Annotation
==========

One of the advantages of the [GWAS-VCF specification](https://github.com/MRCIEU/gwas-vcf-specification) is that summary statistics can be easily annotated.

To obtain a `gff3_file` the following code can be used
```
wget -O- ftp://ftp.ensembl.org/pub/current_gff3/homo_sapiens/Homo_sapiens.GRCh38.108.gff3.gz | gunzip | \
  sed -e 's/^##sequence-region   \([0-9XY]\)/##sequence-region   chr\1/' \
  -e 's/^##sequence-region   MT/##sequence-region   chrM/' \
  -e 's/^\([0-9XY]\)/chr\1/' -e 's/^MT/chrM/' | gzip > $HOME/GRCh38/Homo_sapiens.GRCh38.108.gff3.gz
```

If you want to annotate the coding variants, you can do so with a simple command
```
bcftools csq -Ob -f $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -g $HOME/GRCh38/Homo_sapiens.GRCh38.108.gff3.gz -B 1 -c CSQ -l -n 64 -s - ieu-a-298.hg38.bcf | \
tee ieu-a-298.hg38.csq.bcf | \
bcftools index --force --output ieu-a-298.hg38.csq.bcf.csi
```

You can then quickly extract tables with a list of genome-wide significant variants with coding annotations
```
bcftools +split-vep -Ou -c Consequence -i 'LP>7.3' ieu-a-298.hg38.csq.bcf | \
  bcftools query -f "%CHROM\t%POS[\t%LP]\t%Consequence\n" -i 'LP>7.3'
```

To obtain an `rsid_vcf_file` the following code can be used:
```
wget ftp://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.39.gz{,.tbi}
(echo "##fileformat=VCFv4.2"
bcftools view --header-only GCF_000001405.39.gz | grep ^##INFO=\<ID=RS
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO") > tmp.vcf
(bcftools reheader --fai GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai --temp-prefix ./bcftools. tmp.vcf
bcftools annotate --no-version --remove ID,^INFO/RS GCF_000001405.39.gz | grep -v "^#\|^NT_\|^NW_" | \
sed 's/NC_012920\.[1-9][0-9]*/chrM/;s/^NC_0*\([1-9][0-9]*\)\.[1-9][0-9]*/chr\1/;s/^chr23/chrX/;s/^chr24/chrY/') | \
bcftools norm --no-version --output-type u --multiallelics -any | \
bcftools norm --no-version --output-type u --check-ref w --rm-dup none --fasta-ref GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | \
bcftools sort --output-type b --temp-dir ./bcftools. | \
tee $HOME/GRCh38/GCF_000001405.39.GRCh38.bcf | \
bcftools index --force --output $HOME/GRCh38/GCF_000001405.39.GRCh38.bcf.csi
/bin/rm tmp.vcf
```

Similarly, you can annotate rsID numbers with
```
bcftools annotate --no-version -a $HOME/GRCh38/GCF_000001405.39.GRCh38.bcf -c RS -Ob ieu-a-298.hg38.bcf | \
tee ieu-a-298.hg38.rsid.bcf | \
bcftools index --force --output ieu-a-298.hg38.rsid.bcf
```

Plotting
========

One of the advantages of having summary statistics in a VCF file is the ability to build an index that allows to retrieve and visualize specific regions of interest

Manhattan plot with all available chromosomes
```
assoc_plot.R \
  --cytoband $HOME/GRCh38/cytoBand.txt.gz \
  --vcf GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.hg38.bcf \
  --png GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.png
```
![](GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.png)

If you you generate an annotated version of the summary statistics
```
bcftools csq -Ob -f $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -g $HOME/GRCh38/Homo_sapiens.GRCh38.108.gff3.gz -B 1 -c CSQ -l -n 64 -s - GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.hg38.bcf | \
tee GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.hg38.csq.bcf | \
bcftools index --force --output GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.hg38.csq.bcf.csi
```

You can then plot and highlight in red all variants that are predicted to affect the protein aminoacid sequence
```
assoc_plot.R \
  --cytoband $HOME/GRCh38/cytoBand.txt.gz \
  --vcf GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.hg38.csq.bcf \
  --csq \
  --region chr15 \
  --png GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.chr15.png
```
![](GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.chr15.png)

And you can zoom and plot any region of interest
```
assoc_plot.R \
  --cytoband $HOME/GRCh38/cytoBand.txt.gz \
  --vcf GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.hg38.csq.bcf \
  --csq \
  --region chr15:81413372-86413372 \
  --png GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.adamtsl3.png
```
![](GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.adamtsl3.png)

Examples
========

Attention Deficit Hyperactivity Disorder
----------------------------------------

Download [ADHD summary statistics](https://figshare.com/articles/dataset/adhdSexSpecific2018/19383299) from [2018 ADHD](http://doi.org/10.1016/j.biopsych.2017.11.026) and [2019 ADHD](http://doi.org/10.1038/s41588-018-0269-7) studies
```
wget -O ADHD_female.GCST012597_buildGRCh37.tsv.gz https://figshare.com/ndownloader/files/35310529
wget -O ADHD_male.GCST005362_buildGRCh37.tsv.gz https://figshare.com/ndownloader/files/35310532
wget -O daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta.gz https://figshare.com/ndownloader/files/28169253

for pfx in female.GCST012597 male.GCST005362; do
  zcat ADHD_${pfx}_buildGRCh37.tsv.gz | cut -f1-3,6- |  sed '1 s/orig_//g' | \
  bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s ADHD_${pfx%.GCST0[01][25][35][69][27]}_2018 | \
  bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -Ob | tee ADHD_$pfx.hg38.bcf | \
  bcftools index --force --output ADHD_$pfx.hg38.bcf.csi
done

bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s ADHD_2019 \
  daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta.gz | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -Ob | tee ADHD_2019.hg38.bcf | \
bcftools index --force --output ADHD_2019.hg38.bcf.csi
```

Anxiety Disorder
----------------

Download [ANX summary statistics](https://figshare.com/articles/dataset/panic2019/16602218) from [2019 ANX](http://doi.org/10.1038/s41380-019-0590-2) study
```
wget -O pgc-panic2019.vcf.tsv.gz https://figshare.com/ndownloader/files/30731276

zcat pgc-panic2019.vcf.tsv.gz | sed '/\t$/d' | \
bcftools +munge --no-version -Ou -C colheaders.tsv --fai human_g1k_v37.fasta.fai -s ANX_2019 | \
bcftools +liftover --no-version -Ou -- \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -Ob | tee pgc-panic2019.hg38.bcf | \
bcftools index --force --output pgc-panic2019.hg38.bcf.csi
```

Autism Spectrum Disorder
------------------------

Download [ASD summary statistics](https://figshare.com/articles/dataset/asd2019/14671989) from [2019 ASD](http://doi.org/10.1038/s41588-019-0344-8) study
```
wget -O iPSYCH-PGC_ASD_Nov2017.gz https://figshare.com/ndownloader/files/28169292

bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s ASD_2017 iPSYCH-PGC_ASD_Nov2017.gz | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -Ob | tee ASD_Nov2017.hg38.bcf | \
bcftools index --force --output ASD_Nov2017.hg38.bcf.csi
```

Bipolar Disorder
----------------

Download [BIP summary statistics](https://figshare.com/articles/dataset/PGC3_bipolar_disorder_GWAS_summary_statistics/14102594) from [2021 BIP](http://doi.org/10.1038/s41588-021-00857-4) study
```
wget -O pgc-bip2021-all.vcf.tsv.gz https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/26603681/pgcbip2021all.vcf.tsv.gz
wget -O pgc-bip2021-BDI.vcf.tsv.gz https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/26603690/pgcbip2021BDI.vcf.tsv.gz
wget -O pgc-bip2021-BDII.vcf.tsv.gz https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/26603702/pgcbip2021BDII.vcf.tsv.gz

for pfx in all BDI BDII; do
  zcat pgc-bip2021-$pfx.vcf.tsv.gz | sed '/\t$/d' | \
  bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s BIP_2021_$pfx | \
  bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -Ob | tee pgc-bip2021-$pfx.hg38.bcf | \
  bcftools index --force --output pgc-bip2021-$pfx.hg38.bcf.csi
done
```

Eating Disorders
----------------

Download [ED summary statistics](https://figshare.com/articles/dataset/an2019/14671980) from [2019 ED](http://doi.org/10.1038/s41588-019-0439-2) study
```
wget -O pgcAN2.2019-07.vcf.tsv.gz https://figshare.com/ndownloader/files/28169271

bcftools +munge --no-version -Ou -C colheaders.tsv --fai human_g1k_v37.fasta.fai -s ED_2019 pgcAN2.2019-07.vcf.tsv.gz | \
bcftools +liftover --no-version -Ou -- \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -Ob | tee pgcAN2.2019-07.hg38.bcf | \
bcftools index --force --output pgcAN2.2019-07.hg38.bcf.csi
```

Major Depressive Disorder
-------------------------

Download [MDD summary statistics](https://figshare.com/articles/dataset/mdd2021asi/16989442) from [2021 MDD](http://doi.org/10.1001/jamapsychiatry.2021.2099) study
```
wget -O jamapsy_Giannakopoulou_2021_exclude_whi_23andMe.txt.gz https://figshare.com/ndownloader/files/31424374
wget -O jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb.txt.gz https://figshare.com/ndownloader/files/34437842

for pfx in 23andMe{,_ukb}; do
  bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta \
    -s MDD_2021_$pfx jamapsy_Giannakopoulou_2021_exclude_whi_$pfx.txt.gz | \
  bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -Ob | tee jamapsy_Giannakopoulou_2021_exclude_whi_$pfx.hg38.bcf | \
  bcftools index --force --output jamapsy_Giannakopoulou_2021_exclude_whi_$pfx.hg38.bcf.csi
done
```

OCD & Tourette Syndrome
-----------------------

Download [OCD-TS summary statistics](https://figshare.com/articles/dataset/ts2019/14672232) frm [2019 OCD-TS](http://doi.org/10.1176/appi.ajp.2018.18070857) study
```
wget -O TS_Oct2018.gz https://figshare.com/ndownloader/files/28169940

bcftools +munge --no-version -Ou -C colheaders.tsv --fai human_g1k_v37.fasta.fai -s TS_2018 TS_Oct2018.gz | \
bcftools +liftover --no-version -Ou -- \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -Ob | tee TS_Oct2018.hg38.bcf | \
bcftools index --force --output TS_Oct2018.hg38.bcf.csi
```

Post Traumatic Stress Disorder
------------------------------

Download [PTSD summary statistics](https://figshare.com/articles/dataset/ptsd2019/14672133) from [2019 PTSD](http://doi.org/10.1038/s41467-019-12576-w) study
```
wget -O pts_all_freeze2_overall.results.gz https://figshare.com/ndownloader/files/28169634

bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta \
  -s PTSD_2019 pts_all_freeze2_overall.results.gz | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -Ob | tee pts_all_freeze2_overall.hg38.bcf | \
bcftools index --force --output pts_all_freeze2_overall.hg38.bcf.csi
```

For ancestry specific results on the autosomes
```
wget -O pts_aam_freeze2_overall.results.gz https://figshare.com/ndownloader/files/28169712
wget -O pts_eur_freeze2_overall.results.gz https://figshare.com/ndownloader/files/28169727
wget -O pts_lat_freeze2_overall.results.gz https://figshare.com/ndownloader/files/28169733

echo -e "AFR aam\nEUR eur\nAMR lat" | \
while read anc type; do
  bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta \
    -s PTSD_2019.$anc pts_${type}_freeze2_overall.results.gz | \
  bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -Ob | tee pts_${type}_freeze2_overall.hg38.bcf | \
  bcftools index --force --output pts_${type}_freeze2_overall.hg38.bcf.csi
done
bcftools merge --no-version -m none -Ob pts_{aam,eur,lat}_freeze2_overall.hg38.bcf | \
tee pts_freeze2_overall.hg38.bcf | \
bcftools index --force --output pts_freeze2_overall.hg38.bcf.csi
/bin/rm pts_{aam,eur,lat}_freeze2_overall.hg38.bcf

bcftools +blupx \
  --no-version \
  --beta-cov 1e-7 \
  --samples PTSD_2019.AFR,PTSD_2019.EUR,PTSD_2019.AMR \
  pts_freeze2_overall.hg38.bcf \
  1kg_ldgm.{AFR,EUR,AMR}.bcf \
  --output-type b \
  --log pts_freeze2_overall.hg38.blup1e-7.log | \
tee pts_freeze2_overall.hg38.blup1e-7.bcf | \
bcftools index --force --output pts_freeze2_overall.hg38.blup1e-7.bcf.csi
```

Schizophrenia
-------------

Download [SCZ summary statistics](https://figshare.com/articles/dataset/scz2022/19426775) from [2022 SCZ](http://doi.org/10.1038/s41586-022-04434-5) study
```
wget -O PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz https://figshare.com/ndownloader/files/34517861
wget -O PGC3_SCZ_wave3.primary.chrX.public.v3.vcf.tsv.gz https://figshare.com/ndownloader/files/34517864
wget -O PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv.gz https://figshare.com/ndownloader/files/34517807
wget -O PGC3_SCZ_wave3.core.chrX.public.v3.vcf.tsv.gz https://figshare.com/ndownloader/files/34517825

for type in primary core; do
  for pfx in autosome chrX; do
    zcat PGC3_SCZ_wave3.$type.$pfx.public.v3.vcf.tsv.gz | sed '/\t$/d' | \
    bcftools +munge --no-version -Ou -C colheaders.tsv --fai human_g1k_v37.fasta.fai -s SCZ_2022.$type | \
    bcftools +liftover --no-version -Ou -- \
      -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
    bcftools sort -Ob | tee PGC3_SCZ_wave3.$type.$pfx.public.v3.hg38.bcf | \
    bcftools index --force --output PGC3_SCZ_wave3.$type.$pfx.public.v3.hg38.bcf.csi
  done
  bcftools concat --no-version --allow-overlaps -Ob PGC3_SCZ_wave3.$type.{autosome,chrX}.public.v3.hg38.bcf | \
  tee PGC3_SCZ_wave3.$type.public.v3.hg38.bcf | \
  bcftools index --force --output PGC3_SCZ_wave3.$type.public.v3.hg38.bcf.csi
  /bin/rm PGC3_SCZ_wave3.$type.{autosome,chrX}.public.v3.hg38.bcf{,.csi}
done
```

For ancestry specific results on the autosomes
```
wget -O PGC3_SCZ_wave3.afram.autosome.public.v3.vcf.tsv.gz https://figshare.com/ndownloader/files/34517801
wget -O PGC3_SCZ_wave3.asian.autosome.public.v3.vcf.tsv.gz https://figshare.com/ndownloader/files/34517804
wget -O PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz https://figshare.com/ndownloader/files/34517828
wget -O PGC3_SCZ_wave3.latino.autosome.public.v3.vcf.tsv.gz https://figshare.com/ndownloader/files/34517855

echo -e "AFR afram 9824 5998\nEAS asian 27363 12305\nEUR european 127906 52017\nAMR latino 4324 1234" | \
while read anc type ns nc; do
  bcftools +munge --no-version -Ou -C colheaders.tsv --fai human_g1k_v37.fasta.fai -s SCZ_2022.$anc --ns $ns --nc $nc \
    PGC3_SCZ_wave3.$type.autosome.public.v3.vcf.tsv.gz | \
  bcftools +liftover --no-version -Ou -- \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -Ob | tee PGC3_SCZ_wave3.$type.autosome.public.v3.hg38.bcf | \
  bcftools index --force --output PGC3_SCZ_wave3.$type.autosome.public.v3.hg38.bcf.csi
done
bcftools merge --no-version -m none -Ob PGC3_SCZ_wave3.{afram,asian,european,latino}.autosome.public.v3.hg38.bcf | \
tee PGC3_SCZ_wave3.autosome.public.v3.hg38.bcf | \
bcftools index --force --output PGC3_SCZ_wave3.autosome.public.v3.hg38.bcf.csi
/bin/rm PGC3_SCZ_wave3.{afram,asian,european,latino}.autosome.public.v3.hg38.bcf

bcftools +blupx \
  --no-version \
  --beta-cov 2e-7 \
  --samples SCZ_2022.AFR,SCZ_2022.EAS,SCZ_2022.EUR,SCZ_2022.AMR \
  PGC3_SCZ_wave3.autosome.public.v3.hg38.bcf \
  1kg_ldgm.{AFR,EAS,EUR,AMR}.bcf \
  --output-type b \
  --log PGC3_SCZ_wave3.autosome.public.v3.hg38.blup2e-7.log | \
tee PGC3_SCZ_wave3.autosome.public.v3.hg38.blup2e-7.bcf | \
bcftools index --force --output PGC3_SCZ_wave3.autosome.public.v3.hg38.blup2e-7.bcf.csi
```

Educational Attainment
----------------------

Download [EDU summary statistics](https://thessgac.com/papers/14) from [2022 EDU](http://doi.org/10.1038/s41588-022-01016-z) study
```
wget https://ssgac.s3.amazonaws.com/ReadMe_EA4.txt
wget https://ssgac.s3.amazonaws.com/EA4_additive_p1e-5_clumped.txt
wget https://ssgac.s3.amazonaws.com/EA4_chrX_p1e-5_clumped.txt

for pfx in additive chrX; do
  bcftools +munge --no-version -Ou -C colheaders.tsv --fai human_g1k_v37.fasta.fai -s EA_2022 EA4_${pfx}_p1e-5_clumped.txt | \
  bcftools +liftover --no-version -Ou -- \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -Ob | tee EA4_${pfx}_p1e-5_clumped.hg38.bcf | \
  bcftools index --force --output EA4_${pfx}_p1e-5_clumped.hg38.bcf.csi
done
bcftools concat --no-version --allow-overlaps -Ob EA4_{additive,chrX}_p1e-5_clumped.hg38.bcf | \
tee EA4_p1e-5_clumped.hg38.bcf | \
bcftools index --force --output EA4_p1e-5_clumped.hg38.bcf.csi
/bin/rm EA4_{additive,chrX}_p1e-5_clumped.hg38.bcf{,.csi}
```

Height
------

Download [Height summary statistics](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#2022_GWAS_Summary_Statistics_and_Polygenic_Score_.28PGS.29_Weights) from [2022 Height](http://doi.org/10.1038/s41586-022-05275-y) study
```
wget https://cnsgenomics.com/data/giant_2022/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.gz

bcftools +munge --no-version -Ou -C colheaders.tsv --fai human_g1k_v37.fasta.fai \
  -s HEIGHT_2022 GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.gz | \
bcftools +liftover --no-version -Ou -- \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -Ob | tee GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.hg38.bcf | \
bcftools index --force --output GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.hg38.bcf.csi
```

For ancestry specific results on the autosomes
```
wget https://cnsgenomics.com/data/giant_2022/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_AFR.gz
wget https://cnsgenomics.com/data/giant_2022/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EAS.gz
wget https://cnsgenomics.com/data/giant_2022/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz
wget https://cnsgenomics.com/data/giant_2022/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_HIS.gz
wget https://cnsgenomics.com/data/giant_2022/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_SAS.gz

echo -e "AFR AFR\nEAS EAS\nEUR EUR\nAMR HIS\nSAS SAS" | \
while read anc type; do
  bcftools +munge --no-version -Ou -C colheaders.tsv --fai human_g1k_v37.fasta.fai -s HEIGHT_2022.$anc \
    GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_$type.gz | \
  bcftools +liftover --no-version -Ou -- \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -Ob | tee GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_$type.hg38.bcf | \
  bcftools index --force --output GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_$type.hg38.bcf.csi
done
bcftools merge --no-version -m none -Ob GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_{AFR,EAS,EUR,HIS,SAS}.hg38.bcf | \
tee GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS.hg38.bcf | \
bcftools index --force --output GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS.hg38.bcf.csi
/bin/rm GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_{AFR,EAS,EUR,HIS,SAS}.hg38.bcf

bcftools +blupx \
  --no-version \
  --beta-cov 2e-7 \
  --samples HEIGHT_2022.AFR,HEIGHT_2022.EAS,HEIGHT_2022.EUR,HEIGHT_2022.AMR,HEIGHT_2022.SAS \
  GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS.hg38.bcf \
  1kg_ldgm.{AFR,EAS,EUR,AMR,SAS}.bcf \
  --output-type b \
  --log GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS.hg38.blup2e-7.log | \
tee GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS.hg38.blup2e-7.bcf | \
bcftools index --force --output GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS.hg38.blup2e-7.bcf.csi
```

BMI
---

Download [BMI summary statistics](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#BMI_and_Height_GIANT_and_UK_BioBank_Meta-analysis_Summary_Statistics) from [2018 BMI](http://doi.org/10.1093/hmg/ddy271) study

```
wget https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz

bcftools +munge --no-version -Ou -C colheaders.tsv --fai human_g1k_v37.fasta.fai \
  -s BMI_2018 Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz | \
bcftools +liftover --no-version -Ou -- \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -Ob | tee Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.hg38.bcf | \
bcftools index --force --output Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.hg38.bcf.csi
```

Smoking
-------

Download [Smoking summary statistics](https://conservancy.umn.edu/handle/11299/201564) from [2019 Smoking](http://doi.org/10.1038/s41588-018-0307-5) study

```
wget https://conservancy.umn.edu/bitstream/handle/11299/201564/SmokingInitiation.txt.gz

bcftools +munge --no-version -Ou -C colheaders.tsv --fai human_g1k_v37.fasta.fai -s SMOKING_2019 SmokingInitiation.txt.gz | \
bcftools +liftover --no-version -Ou -- \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -Ob | tee SmokingInitiation.hg38.bcf | \
bcftools index --force --output SmokingInitiation.hg38.bcf.csi
```

Acknowledgements
================

This work is supported by NIH grant [R01 HG006855](http://grantome.com/grant/NIH/R01-HG006855), NIH grant [R01 MH104964](http://grantome.com/grant/NIH/R01-MH104964), NIH grant [R01MH123451](http://grantome.com/grant/NIH/R01-MH123451), and the Stanley Center for Psychiatric Research

Polygenic score examples
========================

<!--ts-->
   * [Examples](#examples)
      * [Attention Deficit Hyperactivity Disorder](#attention-deficit-hyperactivity-disorder)
      * [Anxiety Disorder](#anxiety-disorder)
      * [Autism Spectrum Disorder](#autism-spectrum-disorder)
      * [Bipolar Disorder](#bipolar-disorder)
      * [Anorexia Nervosa](#anorexia-nervosa)
      * [Major Depressive Disorder](#major-depressive-disorder)
      * [Tourette Syndrome](#tourette-syndrome)
      * [Obsessive compulsive disorder](#obsessive-compulsive-disorder)
      * [Hoarding symptoms](#hoarding-symptoms)
      * [Post Traumatic Stress Disorder](#post-traumatic-stress-disorder)
      * [Schizophrenia](#schizophrenia)
      * [Suicide](#suicide)
      * [Educational Attainment](#educational-attainment)
      * [Intelligence](#intelligence)
      * [Height](#height)
      * [BMI](#bmi)
      * [Smoking](#smoking)
      * [Alzheimer](#alzheimer)
      * [Intracranial Volume](#intracranial-volume)
      * [Hippocampal Volume](#hippocampal-volume)
      * [Cortical](#cortical)
      * [Oscillatory Brain Activity](#oscillatory-brain-activity)
      * [Breast Cancer](#breast-cancer)
<!--te-->

Examples
========

Most of the following examples show how to process summary statistics mostly made available through the [Psychiatric Genomics Consortium](http://pgc.unc.edu/for-researchers/download-results/) (PGC) and generate polygenic score loadings

Attention Deficit Hyperactivity Disorder
----------------------------------------

Download [ADHD summary statistics](http://figshare.com/articles/dataset/adhdSexSpecific2018/19383299) [2018 ADHD](http://doi.org/10.1016/j.biopsych.2017.11.026), [2019 ADHD](http://doi.org/10.1038/s41588-018-0269-7) and [2023 ADHD](http://doi.org/10.1038/s41588-022-01285-8) studies
```
wget -O ADHD_female.GCST012597_buildGRCh37.tsv.gz http://figshare.com/ndownloader/files/35310529
wget -O ADHD_male.GCST005362_buildGRCh37.tsv.gz http://figshare.com/ndownloader/files/35310532
wget -O daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta.gz http://figshare.com/ndownloader/files/28169253
wget http://ipsych.dk/fileadmin/iPSYCH/PGC/ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta_2.zip

for pfx in female.GCST012597 male.GCST005362; do
  zcat ADHD_${pfx}_buildGRCh37.tsv.gz | cut -f1-3,6- |  sed '1 s/orig_//g' | \
  bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s ADHD_${pfx%.GCST0[01][25][35][69][27]}_2018 | \
  bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -o ADHD_$pfx.hg38.bcf -Ob --write-index
done

bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s ADHD_2019 \
  daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta.gz | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o ADHD_2019.hg38.bcf -Ob --write-index

unzip -p ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta_2.zip \
  ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.meta.gz | \
bcftools +munge -C colheaders.tsv -f human_g1k_v37.fasta -s ADHD_2023 | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.hg38.bcf -Ob --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 1e-7 \
  --max-alpha-hat2 8e-4 \
  --exclude 'FILTER="IFFY"' \
  ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.hg38.pgs.b1e-7.bcf \ 
  --output-type b \
  --log ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.hg38.pgs.b1e-7.log \
  --write-index
```

Anxiety Disorder
----------------

Download [ANX summary statistics](http://figshare.com/articles/dataset/panic2019/16602218) from [2019 ANX](http://doi.org/10.1038/s41380-019-0590-2) study
```
wget -O pgc-panic2019.vcf.tsv.gz http://figshare.com/ndownloader/files/30731276

zcat pgc-panic2019.vcf.tsv.gz | sed '/\t$/d' | \
bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s ANX_2019 | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o pgc-panic2019.hg38.bcf -Ob --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 2e-7 \
  --max-alpha-hat2 0.008 \
  --exclude 'FILTER="IFFY"' \
  pgc-panic2019.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output pgc-panic2019.hg38.pgs.b2e-7.bcf \
  --output-type b \
  --log pgc-panic2019.hg38.pgs.b2e-7.log \
  --write-index
```

Autism Spectrum Disorder
------------------------

Download [ASD summary statistics](http://figshare.com/articles/dataset/asd2019/14671989) from [2019 ASD](http://doi.org/10.1038/s41588-019-0344-8) study
```
wget -O iPSYCH-PGC_ASD_Nov2017.gz http://figshare.com/ndownloader/files/28169292

bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s ASD_2017 --ns 46351 --nc 18382 iPSYCH-PGC_ASD_Nov2017.gz | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o ASD_Nov2017.hg38.bcf -Ob --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 1e-7 \
  --max-alpha-hat2 8e-4 \
  --exclude 'FILTER="IFFY"' \
  ASD_Nov2017.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output ASD_Nov2017.hg38.pgs.b1e-7.bcf \
  --output-type b \
  --log ASD_Nov2017.hg38.pgs.b1e-7.log \
  --write-index
```

Bipolar Disorder
----------------

Download [BIP summary statistics](http://figshare.com/articles/dataset/bip2024/27216117) from [2024 BIP](http://doi.org/10.1038/s41586-024-08468-9) study
```
wget -O bip2024_afr_no23andMe.gz http://figshare.com/ndownloader/files/49760766
wget -O bip2024_eas_no23andMe.gz http://figshare.com/ndownloader/files/49760769
wget -O bip2024_eur_no23andMe.gz http://figshare.com/ndownloader/files/49760772
# wget -O bip2024_multianc_no23andMe.gz http://figshare.com/ndownloader/files/49760775
# wget -O bip2024_eur_noUKB_no23andMe.gz http://figshare.com/ndownloader/files/49760787
# wget -O bip2024_eur_clinical_no23andMe_update.gz http://figshare.com/ndownloader/files/52065002
# wget -O bip2024_eur_community_no23andMe_noUKB_update.gz http://figshare.com/ndownloader/files/52064996
# wget -O bip2024_eur_community_no23andMe_update.gz http://figshare.com/ndownloader/files/52064999

for pfx in afr eas eur; do
  bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s BIP_2024.${pfx^^} bip2024_${pfx}_no23andMe.gz | \
  bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -o bip2024_${pfx}_no23andMe.hg38.bcf -Ob --write-index
done
bcftools merge --no-version -m none -o bip2024_no23andMe.hg38.bcf -Ob \
  bip2024_{afr,eas,eur}_no23andMe.hg38.bcf --write-index
/bin/rm bip2024_{afr,eas,eur}_no23andMe.hg38.bcf{,.csi}

bcftools +pgs \
  --no-version \
  --beta-cov 5e-8 \
  --max-alpha-hat2 0.001 \
  --samples BIP_2024.AFR,BIP_2024.EAS,BIP_2024.EUR \
  bip2024_no23andMe.hg38.bcf \
  1kg_ldgm.{AFR,EAS,EUR}.bcf \
  --output bip2024_no23andMe.hg38.pgsx.b5e-8.bcf \
  --output-type b \
  --log bip2024_no23andMe.hg38.pgsx.b5e-8.log \
  --write-index
```

Download [BIP summary statistics](http://figshare.com/articles/dataset/PGC3_bipolar_disorder_GWAS_summary_statistics/14102594) from [2021 BIP](http://doi.org/10.1038/s41588-021-00857-4) study
```
wget -O pgc-bip2021-all.vcf.tsv.gz http://s3-eu-west-1.amazonaws.com/pfigshare-u-files/26603681/pgcbip2021all.vcf.tsv.gz
wget -O pgc-bip2021-BDI.vcf.tsv.gz http://s3-eu-west-1.amazonaws.com/pfigshare-u-files/26603690/pgcbip2021BDI.vcf.tsv.gz
wget -O pgc-bip2021-BDII.vcf.tsv.gz http://s3-eu-west-1.amazonaws.com/pfigshare-u-files/26603702/pgcbip2021BDII.vcf.tsv.gz

for pfx in all BDI BDII; do
  zcat pgc-bip2021-$pfx.vcf.tsv.gz | sed '/\t$/d' | \
  bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s BIP_2021_$pfx | \
  bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -o pgc-bip2021-$pfx.hg38.bcf -Ob --write-index

  bcftools +pgs \
    --no-version \
    --beta-cov 2e-7 \
    --max-alpha-hat2 0.001 \
    --exclude 'FILTER="IFFY"' \
    pgc-bip2021-$pfx.hg38.bcf \
    1kg_ldgm.EUR.bcf \
    --output pgc-bip2021-$pfx.hg38.pgs.b2e-7.bcf \
    --output-type b \
    --log pgc-bip2021-$pfx.hg38.pgs.b2e-7.log \
    --write-index
done
```

Anorexia Nervosa
----------------

Download [AN summary statistics](http://figshare.com/articles/dataset/an2019/14671980) from [2019 AN](http://doi.org/10.1038/s41588-019-0439-2) study
```
wget -O pgcAN2.2019-07.vcf.tsv.gz http://figshare.com/ndownloader/files/28169271

bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s AN_2019 pgcAN2.2019-07.vcf.tsv.gz | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o pgcAN2.2019-07.hg38.bcf -Ob --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 2e-7 \
  --max-alpha-hat2 0.001 \
  --exclude 'FILTER="IFFY"' \
  pgcAN2.2019-07.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output pgcAN2.2019-07.hg38.pgs.b2e-7.bcf \
  --output-type b \
  --log pgcAN2.2019-07.hg38.pgs.b2e-7.log \
  --write-index
```

Major Depressive Disorder
-------------------------

Download [MDD summary statistics](http://figshare.com/articles/dataset/GWAS_summary_statistics_for_major_depression_PGC_MDD2025_/27061255) from [2025 MDD](http://doi.org/10.1016/j.cell.2024.12.002) study
```
wget -O pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz http://figshare.com/ndownloader/files/51487019

bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s MDD_2025 pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf -Ob --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 2e-8 \
  --max-alpha-hat2 0.0005 \
  --exclude 'FILTER="IFFY"' \
  pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.pgs.b2e-8.bcf \
  --output-type b \
  --log pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.pgs.b2e-8.log \
  --write-index
```

Download [MDD summary statistics](http://figshare.com/articles/dataset/mdd2021asi/16989442) from [2021 MDD](http://doi.org/10.1001/jamapsychiatry.2021.2099) study (samples sizes estimated from eTable2)
```
wget -O jamapsy_Giannakopoulou_2021_exclude_whi_23andMe.txt.gz http://figshare.com/ndownloader/files/31424374
wget -O jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb.txt.gz http://figshare.com/ndownloader/files/34437842

for pfx in 23andMe{,_ukb}; do
  if [ $pfx == 23andMe ]; then
    ns=98502
    nc=12588
    ne=36886.75
  else
    ns=98003
    nc=12455
    ne=36496.54
  fi
  zcat jamapsy_Giannakopoulou_2021_exclude_whi_$pfx.txt.gz | uniq | \
  sed 's/^\(rs142701510\t8\t2289342\t\)t\tc/\1a\tg/;s/^\(rs566706139\t8\t2289350\t\)t\tg/\1a\tc/;s/^\(rs74664568\t8\t2324610\t\)t\tc/\1a\tg/' | \
  bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta \
    -s MDD_2021 --ns $ns --nc $nc --ne $ne | \
  bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -o jamapsy_Giannakopoulou_2021_exclude_whi_$pfx.hg38.bcf -Ob --write-index
done

bcftools +pgs \
  --no-version \
  --beta-cov 3e-8 \
  --max-alpha-hat2 0.001 \
  --exclude 'FILTER="IFFY"' \
  jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb.hg38.pgs.b3e8.bcf \
  --output-type b \
  --log jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb.hg38.pgs.b3e8.log \
  --write-index
```

Tourette Syndrome
-----------------

Download [TS summary statistics](http://figshare.com/articles/dataset/ts2019/14672232) frm [2019 TS](http://doi.org/10.1176/appi.ajp.2018.18070857) study
```
wget -O TS_Oct2018.gz http://figshare.com/ndownloader/files/28169940

zcat TS_Oct2018.gz | sort -k2,2n -k3,3n | \
sed 's/^\(rs8075185 17 36060216 \)A G/\1T C/;s/^\(rs2855958 7 142170167 \)T C/\1A G/;s/^\(rs17274 7 142224511 \)T C/\1A G/' | \
bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s TS_2019 --ns 14307 --nc 4819 | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o TS_Oct2018.hg38.bcf -Ob --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 2e-7 \
  --max-alpha-hat2 0.002 \
  --exclude 'FILTER="IFFY"' \
  TS_Oct2018.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output TS_Oct2018.hg38.pgs.b2e-7.bcf \
  --output-type b \
  --log TS_Oct2018.hg38.pgs.b2e-7.log \
  --write-index
```

Obsessive compulsive disorder
-----------------------------

Download [OCD summary statistics](http://figshare.com/articles/dataset/ocd2018/14672103) from [2018 OCD](http://doi.org/10.1038/mp.2017.154) study
```
wget -O ocd_aug2017.gz http://figshare.com/ndownloader/files/28169544

bcftools +munge --no-version -Ou -c PLINK -f human_g1k_v37.fasta -s OCD_2018 ocd_aug2017.gz --ns 9725 --nc 2688 | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o ocd_aug2017.hg38.bcf -Ob --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 2e-7 \
  --max-alpha-hat2 0.003 \
  --exclude 'FILTER="IFFY"' \
  ocd_aug2017.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output ocd_aug2017.hg38.pgs.b2e-7.bcf \
  --output-type b \
  --log ocd_aug2017.hg38.pgs.b2e-7.log \
  --write-index
```

Hoarding symptoms
-----------------

Download [hoarding summary statistics](http://figshare.com/articles/dataset/hoarding2022/22177274) from [2022 hoarding](http://doi-org.ezp-prod1.hul.harvard.edu/10.1038/s41398-022-02248-7) study
```
wget -O hoarding2022.vcf.tsv.gz http://figshare.com/ndownloader/files/39399566

zcat hoarding2022.vcf.tsv.gz | grep -v ^$ | \
bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s hoarding_2022 | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o hoarding2022.hg38.bcf -Ob --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 4e-8 \
  --max-alpha-hat2 0.001 \
  --exclude 'FILTER="IFFY"' \
  hoarding2022.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output hoarding2022.hg38.pgs.b4e-8.bcf \
  --output-type b \
  --log hoarding2022.hg38.pgs.b4e-8.log \
  --write-index
```

Post Traumatic Stress Disorder
------------------------------

Download [PTSD summary statistics](http://figshare.com/articles/dataset/ptsd2019/14672133) from [2019 PTSD](http://doi.org/10.1038/s41467-019-12576-w) study
```
wget -O pts_all_freeze2_overall.results.gz http://figshare.com/ndownloader/files/28169634

bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta \
  -s PTSD_2019 pts_all_freeze2_overall.results.gz | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o pts_all_freeze2_overall.hg38.bcf -Ob --write-index
```

For ancestry specific results on the autosomes
```
wget -O pts_aam_freeze2_overall.results.gz http://figshare.com/ndownloader/files/28169712
wget -O pts_eur_freeze2_overall.results.gz http://figshare.com/ndownloader/files/28169727
wget -O pts_lat_freeze2_overall.results.gz http://figshare.com/ndownloader/files/28169733

echo -e "AFR aam\nEUR eur\nAMR lat" | \
while read anc type; do
  bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta \
    -s PTSD_2019.$anc pts_${type}_freeze2_overall.results.gz | \
  bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -o pts_${type}_freeze2_overall.hg38.bcf -Ob --write-index
done
bcftools merge --no-version -m none -o pts_freeze2_overall.hg38.bcf -Ob \
  pts_{aam,eur,lat}_freeze2_overall.hg38.bcf --write-index
/bin/rm pts_{aam,eur,lat}_freeze2_overall.hg38.bcf{,.csi}

bcftools +pgs \
  --no-version \
  --beta-cov 4e-8 \
  --max-alpha-hat2 0.001 \
  --samples PTSD_2019.AFR,PTSD_2019.EUR,PTSD_2019.AMR \
  pts_freeze2_overall.hg38.bcf \
  1kg_ldgm.{AFR,EUR,AMR}.bcf \
  --output pts_freeze2_overall.hg38.pgsx.b4e-8.bcf \
  --output-type b \
  --log pts_freeze2_overall.hg38.pgsx.b4e-8.log \
  --write-index
```

Schizophrenia
-------------

Download [SCZ summary statistics](http://figshare.com/articles/dataset/scz2022/19426775) from [2022 SCZ](http://doi.org/10.1038/s41586-022-04434-5) study
```
wget -O PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz http://figshare.com/ndownloader/files/34517861
wget -O PGC3_SCZ_wave3.primary.chrX.public.v3.vcf.tsv.gz http://figshare.com/ndownloader/files/34517864
wget -O PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv.gz http://figshare.com/ndownloader/files/34517807
wget -O PGC3_SCZ_wave3.core.chrX.public.v3.vcf.tsv.gz http://figshare.com/ndownloader/files/34517825

for type in primary core; do
  for pfx in autosome chrX; do
    zcat PGC3_SCZ_wave3.$type.$pfx.public.v3.vcf.tsv.gz | sed 's/NEFF$/NEFFDIV2/;/\t$/d' | \
    bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s SCZ_2022.$type | \
    bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
      -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
    bcftools sort -o PGC3_SCZ_wave3.$type.$pfx.public.v3.hg38.bcf -Ob --write-index
  done
  bcftools concat --no-version --allow-overlaps -o PGC3_SCZ_wave3.$type.public.v3.hg38.bcf -Ob \
    PGC3_SCZ_wave3.$type.{autosome,chrX}.public.v3.hg38.bcf --write-index
  /bin/rm PGC3_SCZ_wave3.$type.{autosome,chrX}.public.v3.hg38.bcf{,.csi}
done
```

For ancestry specific results on the autosomes
```
wget -O PGC3_SCZ_wave3.afram.autosome.public.v3.vcf.tsv.gz http://figshare.com/ndownloader/files/34517801
wget -O PGC3_SCZ_wave3.asian.autosome.public.v3.vcf.tsv.gz http://figshare.com/ndownloader/files/34517804
wget -O PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz http://figshare.com/ndownloader/files/34517828
wget -O PGC3_SCZ_wave3.latino.autosome.public.v3.vcf.tsv.gz http://figshare.com/ndownloader/files/34517855

echo -e "AFR afram\nEAS asian\nEUR european\nAMR latino" | \
while read anc type; do
  if [ type=="afram" ]; then opt="--ns 9824 --nc 5998 --ne 9234.7"; else opt=""; fi
  if [ type=="latino" ]; then opt="--ns 4324 --nc 1234 --ne 3335.2"; else opt=""; fi
  zcat PGC3_SCZ_wave3.$type.autosome.public.v3.vcf.tsv.gz | sed 's/NEFF$/NEFFDIV2/' | \
  bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s SCZ_2022.$anc $opt | \
  bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -o PGC3_SCZ_wave3.$type.autosome.public.v3.hg38.bcf -Ob --write-index
done
bcftools merge --no-version -m none -o PGC3_SCZ_wave3.autosome.public.v3.hg38.bcf -Ob \
  PGC3_SCZ_wave3.{afram,asian,european,latino}.autosome.public.v3.hg38.bcf --write-index
/bin/rm PGC3_SCZ_wave3.{afram,asian,european,latino}.autosome.public.v3.hg38.bcf{,.csi}

bcftools +pgs \
  --no-version \
  --beta-cov 2e-7 \
  --max-alpha-hat2 0.002 \
  --samples SCZ_2022.AFR,SCZ_2022.EAS,SCZ_2022.EUR,SCZ_2022.AMR \
  --exclude 'FILTER="IFFY"' \
  PGC3_SCZ_wave3.autosome.public.v3.hg38.bcf \
  1kg_ldgm.{AFR,EAS,EUR,AMR}.bcf \
  --output PGC3_SCZ_wave3.autosome.public.v3.hg38.pgsx.b2e-7.bcf \
  --output-type b \
  --log PGC3_SCZ_wave3.autosome.public.v3.hg38.pgsx.b2e-7.log \
  --write-index
```

Suicide
-------

Download [SUI summary statistics](http://tinyurl.com/ISGC2021) from [2022 SUI](http://doi.org/10.1016/j.biopsych.2021.05.029) study (notice that you will need a Dropbox link provided from [PGC.DAC.SUI](mailto:pgc.dac.sui@gmail.com))
```
wget -O daner_model2_062620_eur.neff.qc2.80.gz "http://www.dropbox.com/sh/<ISGC2021_D>/<ISGC2021_ID>_<ISGC2021_ID>/daner_model2_062620_eur.neff.qc2.80.gz?dl=0"

zcat daner_model2_062620_eur.neff.qc2.80.gz | sed 's/\.y\t/\t/g' | \
bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s SUI_2022 | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o daner_model2_062620_eur.neff.qc2.80.bcf -Ob --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 1e-7 \
  --max-alpha-hat2 0.0005 \
  --exclude 'FILTER="IFFY"' \
  daner_model2_062620_eur.neff.qc2.80.bcf \
  1kg_ldgm.EUR.bcf \
  --output daner_model2_062620_eur.neff.qc2.80.pgs.b1e7.bcf \
  --output-type b \
  --log daner_model2_062620_eur.neff.qc2.80.pgs.b1e7.log \
  --write-index
```

Educational Attainment
----------------------

Download [EDU summary statistics](http://thessgac.com/papers/14) from [2022 EDU](http://doi.org/10.1038/s41588-022-01016-z) study
```
wget http://ssgac.s3.amazonaws.com/ReadMe_EA4.txt
wget http://ssgac.s3.amazonaws.com/EA4_additive_p1e-5_clumped.txt
wget http://ssgac.s3.amazonaws.com/EA4_chrX_p1e-5_clumped.txt

for pfx in additive chrX; do
  bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s EA_2022 EA4_${pfx}_p1e-5_clumped.txt | \
  bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -o EA4_${pfx}_p1e-5_clumped.hg38.bcf -Ob --write-index
done
bcftools concat --no-version --allow-overlaps -o EA4_p1e-5_clumped.hg38.bcf -Ob \
  EA4_{additive,chrX}_p1e-5_clumped.hg38.bcf --write-index
/bin/rm EA4_{additive,chrX}_p1e-5_clumped.hg38.bcf{,.csi}
```

Intelligence
------------

Download [IQ summary statistics](http://ctg.cncr.nl/software/summary_statistics/) from [2018 IQ](http://doi.org/10.1038/s41588-018-0152-6) study
```
wget http://ctg.cncr.nl/documents/p1651/SavageJansen_IntMeta_sumstats.zip

unzip -p SavageJansen_IntMeta_sumstats.zip sumstats/SavageJansen_2018_intelligence_metaanalysis.txt | \
bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s IQ_2018 | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o SavageJansen_IntMeta.hg38.bcf -Ob --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 1e-7 \
  --max-alpha-hat2 0.0005 \
  --exclude 'FILTER="IFFY"' \
  SavageJansen_IntMeta.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output SavageJansen_IntMeta.hg38.pgs.b1e7.bcf \
  --output-type b \
  --log SavageJansen_IntMeta.hg38.pgs.b1e7.log \
  --write-index
```

Height
------

Download [Height summary statistics](http://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#2022_GWAS_Summary_Statistics_and_Polygenic_Score_.28PGS.29_Weights) from [2022 Height](http://doi.org/10.1038/s41586-022-05275-y) study
```
wget http://cnsgenomics.com/data/giant_2022/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.gz

bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta \
  -s HEIGHT_2022 GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.gz | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.hg38.bcf -Ob --write-index
```

For ancestry specific results on the autosomes
```
wget http://cnsgenomics.com/data/giant_2022/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_AFR.gz
wget http://cnsgenomics.com/data/giant_2022/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EAS.gz
wget http://cnsgenomics.com/data/giant_2022/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz
wget http://cnsgenomics.com/data/giant_2022/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_HIS.gz
wget http://cnsgenomics.com/data/giant_2022/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_SAS.gz

echo -e "AFR AFR\nEAS EAS\nEUR EUR\nAMR HIS\nSAS SAS" | \
while read anc type; do
  bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s HEIGHT_2022.$anc \
    GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_$type.gz | \
  bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -o GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_$type.hg38.bcf -Ob --write-index
done
bcftools merge --no-version -m none -o GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS.hg38.bcf -Ob \
  GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_{AFR,EAS,EUR,HIS,SAS}.hg38.bcf --write-index
/bin/rm GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_{AFR,EAS,EUR,HIS,SAS}.hg38.bcf{,.csi}

bcftools +pgs \
  --no-version \
  --beta-cov 1e-7 \
  --max-alpha-hat2 0.005 \
  --samples HEIGHT_2022.AFR,HEIGHT_2022.EAS,HEIGHT_2022.EUR,HEIGHT_2022.AMR,HEIGHT_2022.SAS \
  GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS.hg38.bcf \
  1kg_ldgm.{AFR,EAS,EUR,AMR,SAS}.bcf \
  --output GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS.hg38.pgsx.b1e-7.bcf \
  --output-type b \
  --log GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS.hg38.pgsx.b1e-7.log \
  --write-index
```

BMI
---

Download [BMI summary statistics](http://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#BMI_and_Height_GIANT_and_UK_BioBank_Meta-analysis_Summary_Statistics) from [2018 BMI](http://doi.org/10.1093/hmg/ddy271) study

```
wget http://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz

bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta \
  -s BMI_2018 Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.hg38.bcf -Ob --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 2e-7 \
  --max-alpha-hat2 0.002 \
  Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.hg38.pgs.b2e-7.bcf \
  --output-type b \
  --log Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.hg38.pgs.b2e-7.log \
  --write-index
```

Smoking
-------

Download [Smoking summary statistics](http://conservancy.umn.edu/handle/11299/201564) from [2019 Smoking](http://doi.org/10.1038/s41588-018-0307-5) study

```
wget http://conservancy.umn.edu/bitstream/handle/11299/201564/SmokingInitiation.txt.gz

bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s SMOKING_2019 SmokingInitiation.txt.gz | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o SmokingInitiation.hg38.bcf -Ob --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 2e-8 \
  --max-alpha-hat2 0.0005 \
  SmokingInitiation.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output SmokingInitiation.hg38.pgs.b2e-8.bcf \
  --output-type b \
  --log SmokingInitiation.hg38.pgs.b2e-8.log \
  --write-index
```

Alzheimer
---------

Download [Alzheimer summary statistics](http://www.niagads.org/datasets/ng00075) from [2019 Alzheimer](http://doi.org/10.1038/s41588-019-0358-2) study

```
wget -O Kunkle_etal_Stage2_results.txt http://www.niagads.org/system/tdf/public_docs/Kunkle_etal_Stage2_results.txt?file=1

bcftools +munge --no-version -Ou \
  -C colheaders.tsv -f human_g1k_v37.fasta \
  -s AD_2019 Kunkle_etal_Stage2_results.txt | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o Kunkle_etal_Stage2_results.hg38.bcf -Ob --write-index
```

Download [Alzheimer summary statistics](http://www.ebi.ac.uk/gwas/studies/GCST90027158) from [2022 Alzheimer](http://doi.org/10.1038/s41588-022-01024-z) study

```
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/GCST90027158_buildGRCh38.tsv.gz

bcftools +munge --no-version -o GCST90027158.hg38.bcf -Ob -C colheaders.tsv \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -s AD_2022 GCST90027158_buildGRCh38.tsv.gz \
  --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 2e-8 \
  --max-alpha-hat2 0.002 \
  GCST90027158.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --exclude 'FILTER="IFFY"' \
  --output GCST90027158.hg38.pgs.b2e-8.bcf \
  --output-type b \
  --log GCST90027158.hg38.pgs.b2e-8.log \
  --write-index
```

Intracranial Volume
-------------------

Download [ICV summary statistics](http://enigma.ini.usc.edu/research/download-enigma-gwas-results/) from [2016 ICV](http://doi.org/10.1038/nn.4398) study

```
wget http://enigma.ini.usc.edu/wp-content/uploads/E2_C/CHARGE-ENIGMA-ICV-METAANALYSIS-201311141.TBL.FINAL.gz

zcat CHARGE-ENIGMA-ICV-METAANALYSIS-201311141.TBL.FINAL.gz | \
sed 's/^MarkerName/chromosome\tbase_pair_location/;s/:/\t/;/R\t/d;s/ /\t/g' | \
bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s ICV_2016 | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o CHARGE-ENIGMA-ICV-METAANALYSIS-201311141.hg38.bcf -Ob --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 1e-7 \
  --max-alpha-hat2 0.005 \
  CHARGE-ENIGMA-ICV-METAANALYSIS-201311141.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output CHARGE-ENIGMA-ICV-METAANALYSIS-201311141.hg38.pgs.b1e-7.bcf \
  --output-type b \
  --log CHARGE-ENIGMA-ICV-METAANALYSIS-201311141.hg38.pgs.b1e-7.log \
  --write-index
```

Hippocampal Volume
------------------

Download [HV summary statistics](http://enigma.ini.usc.edu/research/download-enigma-gwas-results/) from [2017 HV](http://doi.org/10.1038/ncomms13624) study

```
wget http://enigma.ini.usc.edu/wp-content/uploads/E2_C/CHARGE-ENIGMA-HV-METAANALYSIS-201311141.TBL.FINAL.gz

zcat CHARGE-ENIGMA-HV-METAANALYSIS-201311141.TBL.FINAL.gz | \
sed 's/^MarkerName/chromosome\tbase_pair_location/;s/:/\t/;/R\t/d;s/ /\t/g' | \
bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s HV_2017 | \
bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o CHARGE-ENIGMA-HV-METAANALYSIS-201311141.hg38.bcf -Ob --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 8e-8 \
  --max-alpha-hat2 0.004 \
  CHARGE-ENIGMA-HV-METAANALYSIS-201311141.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output CHARGE-ENIGMA-HV-METAANALYSIS-201311141.hg38.pgs.b8e-8.bcf \
  --output-type b \
  --log CHARGE-ENIGMA-HV-METAANALYSIS-201311141.hg38.pgs.b8e-8.log \
  --write-index
```

Cortical
--------

Download [Cortical summary statistics](http://enigma.ini.usc.edu/downloads/) from [2020 Cereb. Cortex](http://doi.org/10.1126/science.aay6690) study

```
for pfx in SurfArea Thickness; do
  wget http://enigma.ini.usc.edu/downloads/ENIGMA3_Global/ENIGMA3_mixed_se_wo_Mean_Full_${pfx}_20190429.txt.gz

  bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s ${pfx}_2020 \
    ENIGMA3_mixed_se_wo_Mean_Full_${pfx}_20190429.txt.gz | \
  bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -o ENIGMA3_mixed_se_wo_Mean_Full_${pfx}_20190429.hg38.bcf -Ob --write-index
done

bcftools +pgs \
  --no-version \
  --beta-cov 1e-7 \
  --max-alpha-hat2 0.005 \
  ENIGMA3_mixed_se_wo_Mean_Full_SurfArea_20190429.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output ENIGMA3_mixed_se_wo_Mean_Full_SurfArea_20190429.hg38.pgs.b1e-7.bcf \
  --output-type b \
  --log ENIGMA3_mixed_se_wo_Mean_Full_SurfArea_20190429.hg38.pgs.b1e-7.log \
  --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 4e-8 \
  --max-alpha-hat2 0.002 \
  ENIGMA3_mixed_se_wo_Mean_Full_Thickness_20190429.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output ENIGMA3_mixed_se_wo_Mean_Full_Thickness_20190429.hg38.pgs.b4e-8.bcf \
  --output-type b \
  --log ENIGMA3_mixed_se_wo_Mean_Full_Thickness_20190429.hg38.pgs.b4e-8.log \
  --write-index
```

Oscillatory Brain Activity
--------------------------

Download [EEG summary statistics](http://enigma-brain.org/download/sumstats/getFiles.php?key=MmIyZTNmMzljNWRiYWM0OWY2OTAwZGFlNjUxMTNmZTM=
) from [2018 EEG](http://doi.org/10.1002/hbm.24238) study

```
echo alphaCz alphaOcc betaCz deltaCz peakOcc thetaCz | tr ' ' '\n' | cat -n | \
while read n pfx; do
  wget -O ENIGMA-EEG_20181101_$pfx.txt.gz 'http://enigma-brain.org/download/sumstats/getFiles.php?key=MmIyZTNmMzljNWRiYWM0OWY2OTAwZGFlNjUxMTNmZTM=&f='$n

  bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s ${pfx}_2018 \
    ENIGMA-EEG_20181101_$pfx.txt.gz | \
  bcftools +liftover --no-version -Ou -- -s human_g1k_v37.fasta \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
  bcftools sort -o ENIGMA-EEG_20181101_$pfx.hg38.bcf -Ob --write-index
done
```

Breast Cancer
-------------

Download [breast cancer summary statistics](https://drive.google.com/drive/folders/1UZHnFkT8xlcSOPM9ATmpx_QQv2mcQewA?usp=sharing) from [2020 breast cancer](http://doi.org/10.1038/s41588-020-0609-2) study

```
wget -O- CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics.txt https://drive.google.com/file/d/1rMzEjNwyegS3J_eY15szLFfSzi6Vxrwe/view?usp=drive_link

cat CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics.txt | tr -d '"' <  | \
bcftools +munge --no-version -Ou -C colheaders.tsv -f human_g1k_v37.fasta -s BC | \
bcftools +liftover --exclude 'FILTER="REF_MISMATCH"' --no-version -Ou -- -s human_g1k_v37.fasta \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -c hg19ToHg38.over.chain.gz | \
bcftools sort -o CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics.bcf -Ob --write-index

bcftools +pgs \
  --no-version \
  --beta-cov 4e-8 \
  --max-alpha-hat2 0.002 \
  ENIGMA3_mixed_se_wo_Mean_Full_Thickness_20190429.hg38.bcf \
  1kg_ldgm.EUR.bcf \
  --output ENIGMA3_mixed_se_wo_Mean_Full_Thickness_20190429.hg38.pgs.b4e-8.bcf \
  --output-type b \
  --log ENIGMA3_mixed_se_wo_Mean_Full_Thickness_20190429.hg38.pgs.b4e-8.log \
  --write-index
```
# Single Cell Alternative Polyadenylation using Expectation-maximization (SCAPE)

## what is SCAPE?

SCAPE is a tool for estimating alternative polyadenylation (APA) events from single-end or pair-end scRNA-seq data (10x, Microwell-seq). The scRNA-seq must be generated by protocols based on polyT enrichment. SCAPE utilizes the insert size information of DNA fragments to estimate polyadenylation sites (pA sites). 

SCAPE models each pA sites as a Gaussian distribution, the corresponding isoform of which also has a weight (proportion). The final estimations are the mean and standard deviation of all pA sites, as well as their weights. The following figure shows the basic idea of SCAPE.

<img src="https://github.com/zhou-ran/SCAPE/blob/main/image/model_graph.png" alt="model_graph" width="600" />

## Input 

### scRNA-seq data (bam file)
This can be single-end or pair-end sequencing scRNA-seq data. We only tested on data sequenced on Illumina platform. If special treatments have been taken regarding insert size in the library preparation, it is recommended to adjust the insert size distribution in the source code.  

### Gene annotation file (GTF file)
The Gene transfer format (GTF) is a file format used to hold information about gene structure. Since APA is an splicing event in the 3'-UTR, only reads falls in the 3'-UTR is relevant. GTF file is needed to pull out reads that fall into 3'-UTR, which is then used for the analysis.

GTF files for different species can be downloaded from [Ensembl database](https://www.ensembl.org/info/data/ftp/index.html). 

## Output

TODO


## Usage

SCAPE is written in Python. It is a pre-requisite to have python3 installed before running SCAPE. SCAPE has two general functions (1) data preprocessing and (2)  APA event inference. 

```
Usage: python main.py [OPTIONS] COMMAND [ARGS]...

  An analysis framework for analysing alternative polyadenylation at single
  cell levels. Current version: 1.0.0

Options:
  --help  Show this message and exit.

Commands:
  apamix
  prepare

```

### prepare

Preprocessing utr and intron region for inferring alternative polyadenylation events

```
Usage: python main.py prepare [OPTIONS]

Options:
  --gtf TEXT     The gtf (gz or not, but must be sorted) file for preparing utr and intron region.
                 [required]

  --prefix TEXT  The prefix of output bed file.
  --help         Show this message and exit.
```

run

```shell

bedtools sort -i GRCh38.p12.cr.gtf  | bgzip > GRCh38.p12.cr.gtf.gz
python main.py prepare --gtf GRCh38.p12.cr.gtf.gz --prefix GRCh38


```

### apamix

Perform the actual alternative polyadenylation events inference. If one wants to change the default paramters such as insert size distribution, maximum number of pA sites etc, it can be done by modifying the file "apamix.py" under folder "apamix".

```
Usage: python main.py apamix [OPTIONS]

Options:
  --bed TEXT       The target regions (bed format) used for mix inference
                   [required]

  --bam TEXT       The bam file (sorted and indexed)  [required]
  -o, --out TEXT   The output path  [required]
  --cores INTEGER  Num (cores) of region are infering at once
  --cb TEXT        The cell barcode file, one cb for one line.  [required]
  -v, --verbose    Verbose mode
  --help           Show this message and exit.
```

run


```shell

python main.py apamix \
--bed GRCh38_utr.bed \
--bam input.bam \
--out data/ \
--cores 12 \
--cb cellbarcode.tsv

```

## Examples

Visit [wiki](https://github.com/LuChenLab/SCAPE/wiki) to run SCAPE on an example dataset.

## Downstream analysis

TODO
The output of SCAPE can be built into a Seurat object in R, which could be easily use for downstream single cell analyses. 

## Citation

If you use SCAPE in your publication, please cite SCAPE by

Zhou *et al*. SCAPE: A mixture model revealing single-cell polyadenylation diversity and cellular dynamic during cell differentiation and reprogramming.

# twisst2

`twisst2` is a tool for [topology weighting](https://doi.org/10.1534/genetics.116.194720). Topology weighting summarises genealogies in terms of the relative abundance of different sub-tree topologies. It can be used to explore processes like introgression and it can aid the identification of trait-associated loci.

`twisst2` has a number of important improvements over the original [`twisst`](https://github.com/simonhmartin/twisst) tool. Most importantly, `twisst2` **incorporates inference of the ancestral recombination graph (ARG) or tree sequence** - local genealogies and their breakpoints along the chromosome. It does this using [`sticcs`](https://github.com/simonhmartin/sticcs). `sticcs` is a model-free approach and it does not require phased data, so `twisst2` **can run on unphased genotypes of any ploidy**.

The recommended way to run `twisst2` is to start from polarised genotype data. This means you either need to know the ancestral allele at each site, or you need an appropriate outgroup(s) to allow inference of the derived allele.

An alternative way to run it is by first inferring the ancestral recombination graph (ARG) tree sequence using a different tool like [Relate](https://myersgroup.github.io/relate/) or [`tsinfer`](https://tskit.dev/tsinfer/docs/stable/index.html). However, this typically requires phased genotypes, and [my tests](https://doi.org/10.1093/genetics/iyaf181) suggest that `twisst2+sticcs` is more accurate than other methods anyway.

# Publications
- The general concept of topology weighting is described by [Martin and Van Belleghem 2017](https://doi.org/10.1534/genetics.116.194720).
- Combining genealogy inference with `sticcs` and topology weighting with `twisst2` is described by [Martin 2025](https://doi.org/10.1093/genetics/iyaf181).


### Installation

First install [`sticcs`](https://github.com/simonhmartin/sticcs) by following the intructions there.

If you would like to analyse tree sequence objects from tools like [`msprime`](https://tskit.dev/msprime/docs/stable/intro.html) and [tsinfer](https://tskit.dev/tsinfer/docs/stable/index.html), you will also need to install [`tskit`](https://tskit.dev/tskit/docs/stable/introduction.html) yourself. To install `twisst2`:

```bash
git clone https://github.com/simonhmartin/twisst2.git

cd twisst2

pip install -e .
```

### Command line tool

#### Starting from unphased (or phased) genotypes

To perform tree inference and topology weighting, `twisst2` takes as input a modified vcf file that contains a `DC` field, giving the count of derived alleles for each individual at each site.

Once you have a vcf file for your genotype data, make the modified version using `sticcs` (this needs to be installed, see above):
```bash
sticcs prep -i <input vcf> -o <output vcf>  --outgroup <outgroup sample ID>
```

If the vcf file already has the ancestral allele (provided in the `AA` field in the `INFO` section), then you do not need to specifiy outrgoups for polarising.

Now you can run the `twisst2` to count sub-tree topologies:

```bash
twisst2 sticcstack -i <input_vcf> -o <output_prefix> --max_subtrees 512 --ploidy 2 --groups <groupname1> <groupname2> <groupname3> <groupname4> --groups_file
```

#### Starting from pre-inferred trees or ARG (e.g. Relate, tsinfer, argweaver, Singer)

```bash
twisst2 trees -i <input_file> -o <output_prefix> --groups <groupname1> <groupname2> <groupname3> <groupname4> --groups_file
```

### Output

- `<output_prefix>.topocounts.tsv.gz` gives the count of each group tree topology for each interval.
- `<output_prefix>.intervals.tsv.gz` gives the chromosome, start and end position of each interval.

### Full command options

#### Topology weighting from VCF: `twisst2 sticcstack`

 | Command Option | Value type | Description | Required? |
 |-----------------|---------------|------------|-------------|
 | `-i`, `--input_vcf` | FILE | Input VCF file with DC field (make thsi with sticcs prep) | Required |
 | `-o`, `--out_prefix` | CHAR | Output file prefix | Required |
 | `--group_names` | CHAR | Name for each group (separated by spaces) | Required unless `--groups_file` is specified |
 | `--groups` | CHAR | Sample IDs (separated by commas) for each group (separated by spaces) | Required unless `--groups_file` is specified |
 | `--groups_file` | FILE | Optional file with a column for sample ID and group | Required unless `--groups` and `--group_names` are specified |
 | `--max_subtrees` | INT | Maximum number of subtrees to consider; we recommend 512 as a good compromise between speed and accuracy | Required |
 | `--ploidy` | INT | Sample ploidy if all the same. Use `--ploidy_file` if samples differ. | Defaults to 2 if not specified and `--ploidy_file` not specified |
 | `--ploidy_file` | FILE | File with two columns for sample ID and ploidy, separated by whitespace | Required if ploidy differs among individuals |
 | `--output_topos` | FILE | Output file for topologies used | Optional |
 | `--variant_range_only` | | Do not extend first and last trees to the ends of the chromosome; use this for region-specific analyses | Optional |
 | `--unrooted` | | Unroot topologies (results in fewer topologies and therefore fewer counts in the output) | Optional |
 | `--no_second_chances` | | Do not consider SNPs that are separated by incompatible SNPs | Optional, not recommend |
 | `--single_pass` | | Single pass when building trees. Only relevant for ploidy > 1 | Optional, not recommended |
 | `--verbose` | | Verbose output | Optional |

#### Topology weighting from tree sequence: `twisst2 trees`
 
 | Command Option | Value type | Description | Required? |
 |-----------------|---------------|------------|-------------|
 | `-i`, `--input_file` | FILE | Input trees file | Required |
 | `-f`, `--input_format` | CHAR | Input file format {tskit, argweaver, newick} | Required |
 | `-o`, `--out_prefix` | CHAR | Output file prefix | Required |
 | `--chrom_name` | CHAR | Chromosome name for output intervals file | Defaults to 'unknown_chrom' |
 | `--unrooted` | | Unroot topologies (results in fewer topologies and therefore fewer counts in the output) | Optional |
 | `--output_topos` | FILE | Output file for topologies used | Optional |
 | `--group_names` | CHAR | Name for each group (separated by spaces) | Required unless `--groups_file` is specified |
 | `--groups` | CHAR | Sample IDs (separated by commas) for each group (separated by spaces) | Required unless `--groups_file` is specified |
 | `--groups_file` | FILE | Optional file with a column for sample ID and group | Required unless `--groups` and `--group_names` are specified |
 | `--verbose` | | Verbose output | Optional |


### R functions for plotting

Some functions for importing and plotting are provided in the `plot_twisst/plot_twisst.R` script. For examples of how to use these functions, see the `plot_twisst/example_plot.R` script.


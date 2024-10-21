# twisst2

`twisst2` is a tool for [topology weighting](https://doi.org/10.1534/genetics.116.194720). Topology weighting summarises genealogies in terms of the relative abundance of different sub-tree topologies. It can be used to explore processes like introgression and it can aid the identification of trait-associated loci.

`twisst2` has a number of important improvements over the original [`twisst`](https://github.com/simonhmartin/twisst) tool. Most importantly, `twisst2` **incorporates inference of the tree sequence** - local genealogies and their breakpoints along the chromosome. It does this using [`sticcs`](https://github.com/simonhmartin/sticcs). `sticcs` is a model-free approach and it does not require phased data, so `twisst2` **can run on unphased genotypes of any ploidy**.

The standard way to run `twisst2` is to start from polarised genotype data. This means you either need to know the ancestral allele at each site, or you need an appropriate outgroup(s) to allow inference of the derived allele.

An alternative way to run it is by first inferring a tree sequence using a different tool like [`tsinfer`](https://tskit.dev/tsinfer/docs/stable/index.html). However, this requires phased and genotypes, and tests suggests that it is less accurate when sample sizes are small.


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

In standard usage, `twisst2` takes as input a modified vcf file that contains a `DC` field, giving the count of derived alleles for each individual at each site.

Once you have a vcf file for your genotype data, make the modified version using `sticcs` (this needs to be installed, see above):
```bash
sticcs prep -i <input vcf> -o <output vcf>  --outgroup <outgroup sample ID>
```

If the vcf file already has the ancestral allele (provided in the `AA` field in the `INFO` section), then you do not need to specifiy outrgoups for polarising.

Now you can run the `twisst2` to count sub-tree topologies:

```bash
twisst2 -i <input_vcf> -o <output_prefix> --max_iterations 100 --ploidy 2 --groups <groupname1> <groupname2> <groupname3> <groupname4> --groups_file
```

### Output

* `<output_prefix>.topocounts.tsv.gz` gives the count of each group tree topology for each interval.
* `<output_prefix>.intervals.tsv.gz` gives the chromosome, start and end position of each interval.

### R functions for plotting

Some functions for importing and plotting are provided in the `plot_twisst/plot_twisst.R` script. For examples of how to use these functions, see the `plot_twisst/example_plot.R` script.


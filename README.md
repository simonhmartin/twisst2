# twisst2

`twisst2` is a tool for [topology weighting](https://doi.org/10.1534/genetics.116.194720). Topology weighting summarises genealogies in terms of the relative abundance of different sub-tree topologies. It can be used to explore processes like introgression and it can aid the identification of trait-associated loci.

`twisst2` has a number of important improvements over the original [`twisst`](https://github.com/simonhmartin/twisst) tool. Most importantly, `twisst2` **incorporates inference of the tree sequence** - local genealogies and their breakpoints along the chromosome. It does this using [`sticcs`](https://github.com/simonhmartin/sticcs). `sticcs` is a model-free approach and it does not require phased data, so `twisst2` **can run on unphased genotypes of any ploidy**.

The standard way to run `sticcs` is to start from polarised genotype data. This means you either need to know the ancestral allele at each site, or you need an appropriate outgroup(s) to allow inference of the derived allele.

An alternative way to run it is by first inferring a tree sequence using a different tool like [tsinfer](https://tskit.dev/tsinfer/docs/stable/index.html]. However, this requires phased and imputed genotypes, and imputation can introduce biases if model assumptions are violated.


### Installation

First install [`sticcs`](https://github.com/simonhmartin/sticcs)` by following the intructions there.

If you would like to analyse tree sequence objects from tools like [`msprime`](https://tskit.dev/msprime/docs/stable/intro.html) and [tsinfer](https://tskit.dev/tsinfer/docs/stable/index.html], you will also need to install [`tskit`](https://tskit.dev/tskit/docs/stable/introduction.html) yourself..

```bash
git clone https://github.com/simonhmartin/twisst2.git

cd twisst2

pip install -e .
```

### Command line tool

In standard usage, `twisst2` takes as input a modified vcf file that contains a `DC` field, giving the count of derived alleles for each individual at each site.

You can make this using `sticcs`:
```bash
sticcs prep -i <input vcf> -o <output vcf>  --outgroup <outgroup sample ID>
```

If the vcf file already has the ancestral allele (provided in the `AA` field in the `INFO` section), then you do not need to specifiy outrgoups for polarising.

Now you can run the `twisst` to count sub-tree topologies:

```bash
twisst2 -i <input_vcf> -o <output_prefix> --max_iterations 100 --ploidy 2 --groups <groupname1> <groupname2> <groupname3> <groupname4> --groups_file
```
This will make a treesequence file that can be loaded and analysed using [tskit](https://tskit.dev/tskit/docs/stable/introduction.html).



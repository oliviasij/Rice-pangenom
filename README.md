# Rice-pangenome
In this GitHub repository, I present a comprehensive workflow for constructing a pangenome graph tailored to rice genomes. This workflow is designed to capture the complex structural variations and repetitive regions present in rice genomes

Below are a workflow proposition describing the construction of a pangenome for rice.
<img width="670" alt="uml" src="https://github.com/oliviasij/Rice-pangenom/assets/99337399/438e9bbd-2fd9-4e88-9955-a33e841b5b4d">

## Packages
The packages used to construct the pangenome are availble at:
### Minimap2
https://github.com/lh3/minimap2
### Samtools
https://samtools.github.io
### Svim-ASM
https://github.com/eldariont/svim-asm
### bgzip
https://github.com/DataBiosphere/bgzip
### tabix
https://github.com/samtools/tabix
### VG
https://github.com/vgteam/vg
### Truvari
https://github.com/ACEnglish/truvari

### Data
The 12 rice genomes can be found on https://ricerc.sicau.edu.cn/
Short reads can be downloaded at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA820969

## Alignment with minimap
The first step is to align assembly contigs on a reference genome, Nipponbare was used as the reference genome

```minimap2 -a -x asm5 --cs -r2k -t 20 <MSU.fasta> <together.fasta> > alignments.sam```
## Sorting and indexing
Afterwards it is sorted and indexed using samtool
```samtools sort -m4G -@4 -o alignments.sorted.bam alignments.sam```

```samtools index alignments.sorted.bam```
## Svim-ASM
The envelope for svim-asm is activated
```conda activate svimasm_env```
The variant.vcf file is created

```svim-asm haploid svim-asm/ <MSU.fasta> --min_mapq 30 --min_sv_size 50 --max_sv_size 1000```
## Preparation
The file is zipped and indexed

```bgzip -c  svim-asm/variants.vcf > variants.vcf.gz```
```tabix -p vcf variants.vcf.gz```
## Construction of graph
The graph is constructed

```vg construct --threads 64 --progress --handle-sv --alt-paths --reference MSU.fasta --vcf calls.vcf > graph.vg```

Indexing and creating dist file

```vg index --threads 64 --temp-dir ./vgtmp --progress --dist-name graph.dist graph.vg```

Indexing and creating xg file

```vg index --threads 64 --temp-dir ./vgtmp --progress --xg-alts --xg-name graph.xg graph.vg```

Indexing and creating gbwt file

```vg gbwt --num-jobs 64 --temp-dir ./vgtmp --progress --path-cover --xg-name graph.xg --output graph.gbwt```

create gbz file

```vg gbwt --num-jobs 64 --temp-dir ./vgtmp --progress --xg-name graph.xg --graph-name graph.gbz --gbz-format graph.gbwt```

minimizer

```vg minimizer --threads 64 --progress --distance-index graph.dist --output-name graph.min graph.gbz```
## Mapping with giraffe
The short reads are mapped to the pangenome

```vg giraffe --threads 64 --progress --output-format gam --gbz-name graph.gbz --minimizer-name graph.min --dist-name graph.dist --fastq-in  <reads.fastq> > alignments.gam```

creating pack files

```vg pack --threads 64 --xg graph.xg --gam alignments.gam --min-mapq 5 --trim-ends 5 --packs-out alignments.pack```
## Calling variants
using vg call to cal upon variants

```vg call --threads 64 --ploidy 2 --vcf calls.vcf --pack alignments.pack graph.xg > calls-genotyped.vcf```
## Preparation for benchmarking
```bgzip -c svim-asm/variants.vcf > variants.vcf.gz```

```tabix -p vcf variants.vcf.gz```

```bgzip -c calls-genotyped.vcf >calls-genotyped.vcf.gz```

```tabix -p vcf calls-genotyped.vcf.gz```
## Benchmarking
```truvari bench -b variants.vcf.gz -c calls-genotyped.vcf  -o results```

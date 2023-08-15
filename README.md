# Rice-pangenome

minimap2 -a -x asm5 --cs -r2k -t 20 <MSU.fasta> <together.fasta> > alignments.sam

samtools sort -m4G -@4 -o alignments.sorted.bam alignments.sam

samtools index alignments.sorted.bam

conda activate svimasm_env

svim-asm haploid svim-asm/ <MSU.fasta> --min_mapq 30 --min_sv_size 50 --max_sv_size 1000

bgzip -c  svim-asm/variants.vcf > variants.vcf.gz
tabix -p vcf variants.vcf.gz

vg construct --threads 64 --progress --handle-sv --alt-paths --reference MSU.fasta --vcf calls.vcf > graph.vg

vg index --threads 64 --temp-dir ./vgtmp --progress --dist-name graph.dist graph.vg

vg index --threads 64 --temp-dir ./vgtmp --progress --xg-alts --xg-name graph.xg graph.vg

vg gbwt --num-jobs 64 --temp-dir ./vgtmp --progress --path-cover --xg-name graph.xg --output graph.gbwt

vg gbwt --num-jobs 64 --temp-dir ./vgtmp --progress --xg-name graph.xg --graph-name graph.gbz --gbz-format graph.gbwt

vg minimizer --threads 64 --progress --distance-index graph.dist --output-name graph.min graph.gbz

vg giraffe --threads 64 --progress --output-format gam --gbz-name graph.gbz --minimizer-name graph.min --dist-name graph.dist --fastq-in  <reads.fastq> > alignments.gam

vg pack --threads 64 --xg graph.xg --gam alignments.gam --min-mapq 5 --trim-ends 5 --packs-out alignments.pack

vg call --threads 64 --ploidy 2 --vcf calls.vcf --pack alignments.pack graph.xg > calls-genotyped.vcf

bgzip -c svim-asm/variants.vcf > variants.vcf.gz
tabix -p vcf variants.vcf.gz

bgzip -c calls-genotyped.vcf >calls-genotyped.vcf.gz
tabix -p vcf calls-genotyped.vcf.gz

truvari bench -b variants.vcf.gz -c calls-genotyped.vcf  -o results

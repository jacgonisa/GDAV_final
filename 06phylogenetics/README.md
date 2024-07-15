## 1&2.BLAST search
make the database
```
makeblastdb -dbtype prot -in all_reference_proteomes.faa -out 06phylogenetics/all_reference_proteomes.blastdb 2> 06phylogenetics/logs/all_reference_proteomes.err
```
make the search
```
tail -n +2 04de_analysis/results/degs_ids.txt | parallel 'blastp -task blastp -query 05comparative/{}.fa  -db 06phylogenetics/all_reference_proteomes.blastdb  -outfmt 6 -evalue 0.001 > 06phylogenetics/results/{}.blastout'
```

## 3.all sequences with hits (and query sequences)
```
tail -n +2 04de_analysis/results/degs_ids.txt | parallel 'python extract_sequences_from_blast_result.py 06phylogenetics/results/{}.blastout all_reference_proteomes.faa | cat 05comparative/{}.fa - > 06phylogenetics/results/{}_with_homologs.faa'
```
## 4.Build a phylogenetic tree

Alignment
```
tail -n +2 04de_analysis/results/degs_ids.txt | parallel -j1 'mafft 06phylogenetics/results/{}_with_homologs.faa > 06phylogenetics/results/{}_with_homologs.alg 2> 06phylogenetics/logs/{}_with_homologs_mafft.err'
```
Tree

```
tail -n +2 04de_analysis/results/degs_ids.txt | parallel -j1 'iqtree -s 06phylogenetics/results/{}with_homologs.alg -m LG > 06phylogenetics/results/logs/{}.out'
```

too many files in results/, let's move the .treefiles to a new folder

```
mkdir 06phylogenetics/results/trees
cp 06phylogenetics/results/*.treefile 06phylogenetics/results/trees/
```

## 5. Visualize tree

```
python /home/compgenomics/4proteomes/scripts/midpoint_rooting.py 06phylogenetics/results/trees/AQUIFEX_01423_with_homologs.alg.treefile | ete3 view --text
```

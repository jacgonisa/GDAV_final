#FINAL PROJECT#

###INITIAL SET UP

#I am cloning the original environment so I can modify mine to add some functionalities, such as 'parallel'
conda create --name jacob_gdav23 --clone gdav23

mamba install parallel

#SET UP ID

parallel echo {} ::: hightemp normaltemp > ids.txt

# TASK1. Use the tool mOTUs for taxonomic profiling
mkdir 01metagenomics
mkdir 01metagenomics/log
mkdir 01metagenomics/results

cat ids.txt | parallel 'motus profile -f metagenomics-hotspring-{}.1.fq.gz -r metagenomics-hotspring-{}.2.fq.gz -n {} -t 20 -o 01metagenomics/{}.motus 2> 01metagenomics/log/{}.motus.err'

tail -n +4 01metagenomics/hightemp.motus | awk -F"\t" '$2 > 0 {print}' | wc -l

##There are 10 species detected

#Question1.1 Adapted from https://motu-tool.org/tutorial.html:
cat ids.txt | parallel 'tail -n +4 01metagenomics/{}.motus | awk -F"\t" '\''$2 > 0 {print $NF, $0}'\'' | sort -k1,1nr | cut -f 2- -d " " > 01metagenomics/results/{}_report_threshold3.txt'

```
cat 01metagenomics/results/hightemp_report_threshold3.txt 
````
Aquifex aeolicus [ref_mOTU_v31_10705]	0.9050938825
Methanococcus maripaludis [ref_mOTU_v31_01426]	0.0832982117
Desulfofundulus australicus/thermocisternus [ref_mOTU_v31_02213]	0.0017500811
Bacillus vietnamensis [ref_mOTU_v31_04124]	0.0016939191
Bradyrhizobium sp. BTAi1 [ref_mOTU_v31_02667]	0.0016556137
Planctomycetaceae species incertae sedis [meta_mOTU_v31_13845]	0.0016444214
Pseudomonas sp. [ref_mOTU_v31_00131]	0.0016169227
Desulfotomaculum ruminis [ref_mOTU_v31_11808]	0.0014541241
Fusobacterium massiliense [ref_mOTU_v31_10166]	0.0014541241
unassigned	0.0003386997

##Aquifex aeolicus is the most abundant, whose taxon id is 224324

# Question1.2
##I am just practising to get info directly from terminal, but a web-based search was done anyways
efetch -db taxonomy -id 224324 -format xml

# Question1.3
##I am just practising to get info directly from terminal, but a web-based search was done anyways
efetch -db nuccore -id NC_000918.1 -format gb > 01metagenomics/NC_000918.1.gb

# Question1.4
cat ids.txt | parallel 'motus profile -f metagenomics-hotspring-{}.1.fq.gz -r metagenomics-hotspring-{}.2.fq.gz -n {} -g 4 -t 20 -o 01metagenomics/{}_threshold4.motus 2> 01metagenomics/logs/{}_threshold4.motus.err'

```
cat ids.txt | parallel 'tail -n +4 01metagenomics/{}_threshold4.motus | awk -F"\t" '\''$2 > 0 {print $NF, $0}'\'' | sort -k1,1nr | cut -f 2- -d " " > 01metagenomics/results/{}_report_threshold4.txt'

cat 01metagenomics/results/hightemp_report_threshold4.txt 

```

Aquifex aeolicus [ref_mOTU_v31_10705]	0.9154098245
Methanococcus maripaludis [ref_mOTU_v31_01426]	0.0842476155
unassigned	0.0003425600

## Still the same, also when threshold=5

cat 01metagenomics/results/hightemp_report_threshold4.txt 
Aquifex aeolicus [ref_mOTU_v31_10705]	0.9154098245
Methanococcus maripaludis [ref_mOTU_v31_01426]	0.0842476155
unassigned	0.0003425600



# Question1.5
##The table was already generated in question1.1 indeed
cat 01metagenomics/results/normaltemp_report_threshold3.txt 
quifex aeolicus [ref_mOTU_v31_10705]	0.0333818811
Methanococcus maripaludis [ref_mOTU_v31_01426]	0.0314499527
Streptococcus pneumoniae [ref_mOTU_v31_00282]	0.0091936821
Clostridiales species incertae sedis [meta_mOTU_v31_12841]	0.0055220382
Oscillibacter sp. 57_20 [meta_mOTU_v31_12610]	0.0054808489
Desulfallas alcoholivorax [ref_mOTU_v31_06842]	0.0054004846
...
## Again, Aquifex aeolicus


# Question1.7
##Just counting the number of species in normal temperature

tail -n +4 01metagenomics/normaltemp.motus | awk -F"\t" '$2 > 0 {print}' | wc -l

##191

# Question1.8
```
paste <(cat normaltemp_report_threshold3.txt | awk '{print $1}' | parallel -j1 'term={}; result=$(esearch -db taxonomy -query "$term" -tool edirect | efetch -format xml -mode xml | xmllint --xpath "/TaxaSet/Taxon/LineageEx/Taxon[Rank/text() = '\''phylum'\'']/ScientificName/text()" - && echo '' || echo "ERROR"); echo -e "$result"') normaltemp_report_threshold3.txt > normaltemp_with_phylum.txt
```

###This command prints the phylum taking the species name as input. No algae eukaryotes were detected (some eukaryotes are detected but I realized this is because a genus name was similar to an insect's, so esearch mislabeled it).
###Further evidence in the report (.pdf)





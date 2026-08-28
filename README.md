Building KrakenUniq database
===============================================

This repository builds a custom KrakenUniq reference database at k-mer length 22 for ancient DNA taxonomic classification. It uses the RefSeq release-232 genomes and follows the Vernot et al. 2021 recipe.

#### **1. Why build a new database**
---------------------
The published and prebuilt databases were built using a k-mer size of 35. Ancient reads are short and deaminated, so a long exact k-mer often finds no match and the majority of the reads are unclassified. A shorter k (~22-25) have potential to recover them. 
KrakenUniq reports distinct k-mers + coverage per taxon, which gives a good control to evaluate real low-abundance signal from a false positive.


- KrakenUniq 1.0.4
- k-mer size = 22
- Minimizer = 15
- DUST-mask + hard-mask to N
- Taxids injected into FASTA headers
- Viral genomes not used

#### **2. Input**
-----------------------------------------------------------------------
Produced by steps 1-2 of prepare_inputs.txt
| File name | Description|
|----------- | -----------|
|assembly_summary_<group>.txt | NCBI tables. col1 = accession, col6 = TAXID, col20 = ftp path. |
|ftp_paths_*.txt|Download URLs (one per assembly)|
|fastas/*_genomic.fna.gz| Downloaded genomes|

#### **3. Build scripts**
--------------------

| File name | Description|
|----------- | -----------|
|build_krakenuniq_db.sbatch | SLURM wrapper: Reserves the node, fixes paths so the KrakenUniq helper scripts are found, sets threads, and runs the .sh from the submit directory.|
|build_krakenuniq_db.sh| the recipe (steps 0-5 below). Reads the inputs, builds the DB, optionally copies it out. Prepared for running a test set, and a full run. |

How to run:
```bash
mkdir -p ~/logs
cd <this folder>

# Test run first (~50 genomes included):
PILOT_N=50 DB=$PWD/DB_pilot sbatch --time=05:00:00 build_krakenuniq_db.sbatch

# Full run (after the pilot is OK)
DELETE_GZ=1 PERSIST_DIR=/ADD/YOUR/STORAGE/PATH/krakenuniq_refseq232_k22 sbatch build_krakenuniq_db.sbatch # with PERSIST_DIR copies result to permanent storage (optional)
```
Parameters defined:
- PILOT_N
- KMER (22)
- THREADS
- NSHARDS
- DELETE_GZ
- PERSIST_DIR
- DB
- SCR

#### **4. Steps descriptions and outputs**
--------------------------------------

| Steps | Description| Message printed in .out log |
|----------- | -----------| -----------|
|0. Check inputs|confirm fastas/*.gz + assembly_summary_*.txt exist and tools on path|  "inputs OK: N genomes, M summaries"
|1. Taxonomy | krakenuniq-build --download-taxonomy (skipped if already present)| "Download taxdump.tar.gz" (saved under DB/taxonomy/)
|2. Taxid map|read the summaries into an accession - taxid table| "taxid map: N accessions"
|3. Mask + tag + add | Per genome: dustmasker masks low-complexity; awk turns masked bases to N (hard-mask) and adds >kraken:taxid <taxid> into headers. Genomes are concatenated into NSHARDS shards; each shard is added to the library and then deleted. (skipped if /library/added already has files) | "masked n/N"; "add-to-library shard_N.fna" (genomes saved in DB/library/added/)
|4. Build index | krakenuniq-build --build --kmer-len 22 (high RAM) | database.kdb, database.idx, database.kdb.counts, taxDB added to DB/
|5. Use STORAGE_DIR (optional)|rsync the runtime DB (excluding library/) to PERSIST_DIR | Add DB to the storage folder

> [!CAUTION] 
> Only use this if the storage folder can be used as input foler for running the pipeline - be aware of network pressure!!

What --build does internally (KrakenUniq stages):
```
database.jdb (jellyfish k-mer count)
|
db_sort
|
database0.kdb
|
seqid2taxid.map 
|
taxDB
|
build LCA database
```

Final verifications that can be done:
```bash
$ krakenuniq-report  --db <DB> | grep -E 'Homo sapiens|Mammalia'
$ krakenuniq-inspect --db <DB> | head
# Seeing 9606 Homo sapiens + mammalian taxids confirms the taxid addition worked.
```
---------------------------------------
Notes
---------------------------------------
- KrakenUniq 1.0.4's krakenuniq-build calls helper scripts via Perl exec. The conda package puts them under share/.../libexec/, not under bin/, so the .sbatch adds that libexec dir to path. Without it, fails.
- krakenuniq-build --add-to-library can exits 255. The .sh tolerates it and checks the file landed in DB/library/added/.
- The build uses --jellyfish-hash-size 4G --work-on-disk to keep the k-mer counting memory limited. Raise the hash size accordingly if the pilot has room for RAM.
- .sh prints the failing line + path so failures are not silent (with ERR trap).


Usage of built DB in sedaDNA_pipeline:
-------------------------------------------
```
workflow/envs/classifier.yaml -> uncomment - krakenuniq=1.0.4

config/config.yaml:
    taxonomy.tool: krakenuniq
    taxonomy.databases.krakenuniq: path/to/database/directory
    paths.results_dir: results_krakenuniq_refseq232_k22
```

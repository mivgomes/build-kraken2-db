Kraken2 custom database
===============================================

This repository builds a custom Kraken2 reference database at k-mer length 22 for ancient DNA taxonomic classification. It uses the RefSeq release-232 genomes and follows the Vernot et al. 2021 study.

#### **1. Why build a new database**
---------------------
The [published and prebuilt databases](https://benlangmead.github.io/aws-indexes/k2) were built using a k-mer with length 35. Ancient reads are short and deaminated, so a long k-mer can find no match and the majority of the reads are unclassified. A shorter k (~22-25) has the potential to recover them.
Kraken2 does not report distinct k-mers per taxon like KrakenUniq, but [`--report-minimizer-data`](https://github.com/DerrickWood/kraken2/wiki/Manual/583de2058f31275b98e3865d0e5a43f906b650dd#distinct-minimizer-count-information) gives distinct-minimizer counts, and [`--confidence`](https://github.com/DerrickWood/kraken2/wiki/Manual/583de2058f31275b98e3865d0e5a43f906b650dd#confidence-scoring) 
adjusts the level of classification, and these two parameters together can serve as a false-positive control at classification time.

- Kraken2
- k-mer size (k) = 22 (follwed k-mer size used in [quicksand](https://olivierrue.pages-forge.inrae.fr/mypage/posts/kraken2_reduce_kmer_length/))
- Minimizer length (l) = 22   (l = k, so every k-mer is its own minimizer, followed approach used by [Rué O. 2026](https://github.com/mpieva/quicksand-build))
- Minimizer spaces (s) = 0 (Default is 7, and if s>0 it cuts the match to l-s informative bases)
- Sequences are masked for low complexity using dustmasker (kraken2 build runs with `--no-masking`, because pre-masking was already done)
- Taxids added into FASTA headers
- Viral genomes not used

#### **2. Input**
-----------------------------------------------------------------------
RefSeq release 232 pre downloaded. Latest release can be found [here](https://ftp.ncbi.nlm.nih.gov/refseq/release/).

| File name | Description|
|----------- | -----------|
|assembly_summary_<group>.txt | NCBI tables. col1 = accession, col6 = TAXID, col20 = ftp path. |
|ftp_paths_*.txt|Download URLs (one per assembly)|
|fastas/*_genomic.fna.gz| Downloaded genomes|

#### **3. Build scripts**
--------------------

| File name | Description|
|----------- | -----------|
|build_kraken2_db.sbatch | SLURM wrapper: reserves the node, puts kraken2-build / dustmasker / GNU parallel on PATH, sets threads, and runs the .sh from the submit directory.|
|build_kraken2_db.sh| the recipe (steps 0-5 below). Reads the inputs, builds the DB, optionally copies it out. Prepared for a test set and a full run. |
|mask_one.sh| helper called by GNU parallel: masks + hard-masks + taxid-tags **one** genome into a shard. A single genome that fails is logged and skipped (exit 0), so it never aborts the whole parallel run. |

How to run:
```bash
mkdir -p ~/logs
cd <this folder>

# Test run first (~50 genomes incl. human GRCh38 + mouse, throwaway DB):
TEST_N=50 DB=$PWD/DB_test sbatch build_kraken2_db.sbatch

# Full run (after the test is OK) - use a FRESH empty DB:
DELETE_GZ=1 STORAGE_DIR=/ADD/YOUR/STORAGE/PATH/kraken2_refseq232_k22 sbatch build_kraken2_db.sbatch
```
Parameters (default):
- TEST_N (0), builds only the first N genomes
- KMER (22)
- MINIMIZER (22)
- THREADS (32), parallel masking jobs + build threads. It comes from the slurm `--cpus-per-task`.
- ROUND_SIZE (2000), genomes masked per round before add-to-library (bounds the local space peak)
- DELETE_GZ (0), add this parameter only if need local space - it deletes the fna.gz files after masking (optional)
- STORAGE_DIR, rsync the finished database in the selected path (optional)
- DB, where the database is built

> [!TIP]
> `DB` and `STORAGE_DIR` are env vars, so a test and a full run can coexist with distinct paths without editing the script. The full run must always target a **fresh empty DB** - see the note on `SKIP_LIBRARY` below.

#### **4. Steps descriptions and outputs**
--------------------------------------

| Steps | Description| Message printed in .out log |
|----------- | -----------| -----------|
|0. Check inputs|confirm fastas/*.gz + assembly_summary_*.txt exist and kraken2-build / dustmasker / parallel + mask_one.sh are present| "inputs OK: N genomes, M summaries"
|1. Taxonomy | `kraken2-build --download-taxonomy --use-ftp --skip-maps` (skipped if already present). `--skip-maps` avoids downloading the accession2taxid maps, taxids already added to the headers.| taxonomy saved under DB/taxonomy/ (names.dmp, nodes.dmp)
|2. Taxid map|read the summaries into an accession -> taxid table| "taxid map: N accessions"
|3. Mask + tag + add (**parallel**) | Per genome: dustmasker masks low-complexity; awk turns masked bases to N (hard-mask) and injects `>kraken:taxid\|<taxid>\|` into headers (via `mask_one.sh`). Masking runs with parallel over per-slot shards, in ROUND_SIZE-genome rounds; each round's shards are then added to the library (`--no-masking`) and deleted. (whole step skipped if /library/added already has files) | "round R: masking N genomes"; "add-to-library shard_S.fna" (genomes saved in DB/library/added/)
|4. Build index | `kraken2-build --build --kmer-len 22 --minimizer-len 22 --minimizer-spaces 0 --no-masking --max-db-size 1000000000000` (minimizer-based; high RAM peak) | hash.k2d, opts.k2d, taxo.k2d added to DB/
|5. Use STORAGE_DIR (optional)|rsync the runtime DB (excluding library/) to STORAGE_DIR | Add DB to the storage folder

> [!CAUTION]
> Only use STORAGE_DIR if the storage folder can be used as input folder for running the pipeline - be aware of network pressure!!

Final verifications that can be done:
```bash
$ kraken2-inspect --db <DB> | grep -E 'Homo sapiens|Mammalia'
# Seeing 9606 Homo sapiens + mammalian taxids confirms the taxid injection worked.
# There should be NO *.fna.masked files under DB/library/ (--no-masking must be on --add-to-library, so kraken2 does not dustmask a second time).
```

Building a Kraken2 database
===============================================

This repository builds a custom Kraken2 reference database at k-mer length 22 for ancient DNA taxonomic classification. It uses the RefSeq release-232 genomes and follows the Vernot et al. 2021 recipe. It is the **Kraken2 twin** of the KrakenUniq k=22 build (`build-krakenuniq-db/`), so the benchmark can separate the *k-mer* effect from the *tool* effect on the same genomes.

#### **1. Why build a new database**
---------------------
The published and prebuilt databases (e.g. `k2_core_nt`, k=35) were built using a long k-mer. Ancient reads are short and deaminated, so a long exact k-mer often finds no match and the majority of the reads are unclassified. A shorter k (~22-25) has the potential to recover them.
Kraken2 does not report distinct k-mers per taxon like KrakenUniq, but `--report-minimizer-data` gives distinct-minimizer counts that serve a similar false-positive control at classification time.

- Kraken2
- k-mer size (k) = 22
- Minimizer length (l) = 20
- Minimizer spaces (s) = 5   (spaced seeds tolerate mismatches, which helps C→T-damaged reads)
- DUST-mask + hard-mask to N (pre-masked, so the build runs with `--no-masking`)
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
|build_kraken2_db.sbatch | SLURM wrapper: reserves the node, puts kraken2-build / dustmasker / GNU parallel on PATH, sets threads, and runs the .sh from the submit directory.|
|build_kraken2_db.sh| the recipe (steps 0-5 below). Reads the inputs, builds the DB, optionally copies it out. Prepared for a test set and a full run. |
|mask_one.sh| helper called by GNU parallel: masks + hard-masks + taxid-tags **one** genome into a shard. A single genome that fails is logged and skipped (exit 0), so it never aborts the whole parallel run. |

How to run:
```bash
mkdir -p ~/logs
cd <this folder>

# Test run first (~50 genomes incl. human GRCh38 + mouse, throwaway DB):
PILOT_N=50 DB=$PWD/DB_pilot sbatch build_kraken2_db.sbatch

# Full run (after the pilot is OK) — use a FRESH empty DB:
DELETE_GZ=1 STORAGE_DIR=/ADD/YOUR/STORAGE/PATH/kraken2_refseq232_k22 sbatch build_kraken2_db.sbatch
```
Parameters (all env-driven, with defaults):
- PILOT_N (0) — >0 builds only the first N genomes
- KMER (22)
- MINIMIZER (20) — l, must be ≤ k and ≤ 31
- SPACES (5) — s, small relative to l
- THREADS — parallel masking jobs + build threads (from `SLURM_CPUS_PER_TASK`)
- ROUND_SIZE (2000) — genomes masked per round before add-to-library (bounds the scratch peak)
- DELETE_GZ (0) — delete staged *.gz after masking
- STORAGE_DIR — rsync the finished runtime DB here (optional)
- DB — where the DB is built (`SHARD_DIR` derives from it as `$DB.shards`)
- SCR — where the staged inputs live

> [!TIP]
> `DB` and `STORAGE_DIR` are env vars, so a pilot and a full run can coexist with distinct paths without editing the script (`SHARD_DIR` follows `DB` automatically). The full run must always target a **fresh empty DB** — see the note on `SKIP_LIBRARY` below.

#### **4. Steps descriptions and outputs**
--------------------------------------

| Steps | Description| Message printed in .out log |
|----------- | -----------| -----------|
|0. Check inputs|confirm fastas/*.gz + assembly_summary_*.txt exist and kraken2-build / dustmasker / parallel + mask_one.sh are present| "inputs OK: N genomes, M summaries"
|1. Taxonomy | `kraken2-build --download-taxonomy --use-ftp --skip-maps` (skipped if already present). `--skip-maps` avoids the multi-GB accession2taxid maps — taxids come from the injected headers, not those maps.| taxonomy saved under DB/taxonomy/ (names.dmp, nodes.dmp)
|2. Taxid map|read the summaries into an accession → taxid table| "taxid map: N accessions"
|3. Mask + tag + add (**parallel**) | Per genome: dustmasker masks low-complexity; awk turns masked bases to N (hard-mask) and injects `>kraken:taxid\|<taxid>\|` into headers (via `mask_one.sh`). Masking runs **THREADS-way in parallel** over per-slot shards, in ROUND_SIZE-genome rounds; each round's shards are then added to the library **serially** (`--no-masking`) and deleted. (whole step skipped if /library/added already has files) | "round R: masking N genomes"; "add-to-library shard_S.fna" (genomes saved in DB/library/added/)
|4. Build index | `kraken2-build --build --kmer-len 22 --minimizer-len 20 --minimizer-spaces 5 --no-masking` (minimizer-based; high RAM peak) | hash.k2d, opts.k2d, taxo.k2d added to DB/
|5. Use STORAGE_DIR (optional)|rsync the runtime DB (excluding library/) to STORAGE_DIR | Add DB to the storage folder

> [!CAUTION]
> Only use STORAGE_DIR if the storage folder can be used as input folder for running the pipeline — be aware of network pressure!!

Final verifications that can be done:
```bash
$ kraken2-inspect --db <DB> | grep -E 'Homo sapiens|Mammalia'
# Seeing 9606 Homo sapiens + mammalian taxids confirms the taxid injection worked.
# There should be NO *.fna.masked files under DB/library/ (--no-masking must be on
# --add-to-library, so kraken2 does not dustmask a second time).
```
---------------------------------------
Notes
---------------------------------------
- Kraken2 is self-contained: **no** libexec PATH hack and **no** exit-255 handling (both were needed for KrakenUniq). The only extra PATH requirement is **GNU parallel** — it must be on PATH for the parallel masking step (add its conda env bin in the .sbatch if it lives elsewhere).
- Masking (step 3) is the slow bottleneck (one dustmasker per genome). Running it THREADS-way in parallel is the main speedup over the serial KrakenUniq build; the index build (step 4) is lighter/faster for Kraken2.
- add-to-library stays **serial** on purpose: concurrent `--add-to-library` races on shared DB state (`library/added`, `prelim_map.txt`). Only the masking is parallelized.
- `SKIP_LIBRARY` is all-or-nothing: if `DB/library/added` has *any* files, the whole masking+add step is skipped and only the index is (re)built. This lets you re-run the build after an OOM at step 4 without redoing masking — but it also means a **full run must never point at a partial DB**, or it would build an incomplete index. Always use a fresh empty DB for a full run.
- .sh prints the failing line + PATH so failures are not silent (ERR trap).

Usage of built DB in sedaDNA_pipeline:
-------------------------------------------
```
workflow/envs/classifier.yaml -> uncomment - kraken2

config/config.yaml:
    taxonomy.tool: kraken2
    taxonomy.databases.kraken2: path/to/database/directory
    paths.results_dir: results_kraken2_refseq232_k22
```

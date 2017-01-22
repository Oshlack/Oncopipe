This directory contains each batch.

From the oncopipe base directory, create the batch and data directory:
```
mkdir -p batches/{BATCH_NAME}/data
```

Place data files (or symlink) in {BATCH_NAME}/data
```
cp /path/to/data/file {BATCH_NAME}/data
```

From the oncopipe base directory, create the batch, and sample.txt for analysis
```
./pipeline/scripts/create_batch.sh {BATCH_NAME} B_ALL designs/B_ALL/B_ALL.bed
```

Run analysis from {BATCH_NAME}/analysis
```
cd {BATCH_NAME}/analysis
bpipe run ../../../pipeline/pipeline.groovy ../samples.txt
```

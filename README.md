# Brcaness-classifier

### To setup Docker
Use documentation [here](https://docs.docker.com/get-docker/) or provided by your Linux distribution.

## Building the docker image
##### dependent on your docker configuration, you may have to add the flag '--network=host'
##### Takes a few hours, is required once.

```bash
tag=0.7
docker build --tag=ovabrca:$tag .
```

======

## Data access for the container

One should never link sensitive parts of a filesystem in a docker container. Read about security [here](https://docs.docker.com/engine/security/security/#docker-daemon-attack-surface).

### A configuration and empty output directory are required
```bash
mkdir -p config/ output/
```

### If the fastq files are in the current directory, then use:

```bash
input="-v `pwd`:/input"
```

### A docker container cannot follow symlinks. You can hard link files - if they are on the same physical disk.
### Otherwise you can (bind) mount the files. e.g. in an input/ directory

```bash
mkdir -p input/
input="-v `pwd`/input:/input"
```

### Even if fastq files were sequenced under the same filename across runs, they can be linked / mounted via input subdirectories.

### to add fastq files, symlinked in current directory, with raw data in a standard Illumina directory structure:

```bash
while read d; do
  p=input/$(basename $d)
  mkdir $p
  input="$input -v $d:/$p:ro"
done < <(ls -1 *.fastq.gz | xargs -n 1 readlink | xargs -n 1 dirname | xargs -n 1 dirname | sort -u)
```

## Configuration files

### To resolve file(s) per sample a tab-separated file containing regular expression and sample basename:
**config/files.txt**

### example
```
.*5610_89_CF25885_.*.gz 5610_89_CF25885
.*5610_8_CF13235_.*.gz  5610_8_CF13235
```

### For the pipeline settings:
**config/config.txt**

### example
```
MEM=3G
THREADS=4
PARALLEL=2
TYPE=2
BRCA_NUM=1
BLACKLIST=GRCh38-blacklist-merged.bed
TAG=-GRCh38-blacklist-merged
```
#### TYPE and BRCA_NUM are currently (version 0.7) unused

# To continue e.g. when more files were added later, you can use:
```
CONTINUE=true
CHECKSUM=false
```

======

## Running the container
### path to bwa indexed fasta files

bwaindex_dir=$path_to/Homo_sapiens.GRCh38/index/bwa/

# full paths are required in docker
path="`pwd`"

### This runs the container 
docker run --rm -u $UID:$GROUPS \
  -v $path/output/:/output:rw \
  $input \
  -v $bwaindex_dir:/ref:ro \
  -v $path/config:/config:ro \
  ovabrca:$tag

## container development
### If possible without docker rebuilding, one can mount the app/ directory, adapt and run scripts manually.
docker run --rm -t -i -u $UID:$GROUPS \
  -v $path/output/:/output:rw \
  $input \
  -v $bwaindex_dir:/ref:ro \
  -v $path/config:/config:ro \
  -v $path/:/app:rw \
  ovabrca:$tag bash


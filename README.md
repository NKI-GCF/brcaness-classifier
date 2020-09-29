# Brcaness-classifier

To setup Docker use documentation [here](https://docs.docker.com/get-docker/) or provided by your Linux distribution.

## Building the docker image
Dependent on your docker configuration, you may have to add the flag '--network=host'.

```bash
tag=0.9
docker build --tag=ovabrca:$tag .
```

======
## Preparing docker container requirements

When the docker image is run, the -v option allows to bind mount a directory (or a file) to a docker container. This docker requires a ref and output mounts and input, unless continuing with previous data.

Docker does not follow file permissions, so do not expose sensitive parts of a filesystem to a docker container. Read about security [here](https://docs.docker.com/engine/security/security/#docker-daemon-attack-surface).

The ref mount should contain a fasta that is bwa indexed, meaning initialized for bwa. The docker image will index, if the fasta was not already.


```bash
cd $your_dedicated_genome_directory
wget ftp://ftp.ensembl.org/pub/release-76/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primaryassembly.fa.gz
reference="-v `pwd`:/ref:rw"
cd -
```

## Data access for the container

An input bind mount is required if processing from raw data.

```bash
mkdir -p output/ input/
input="-v $path_to_fastq_files:/input:ro"
```

If intermediate files (cram files, counts or the classifier text files) are present and survive consistancy checks, processing is continued from there.

A docker container cannot follow symlinks. You can hard link files - if they are on the same physical disk. For samples sequenced with the same name in different runs docker mounts can be nested in input/ subdirectories; e.g.:

```bash
while read d; do
  p=input/$(basename $d)
  mkdir $p
  input="$input -v $d:/$p:ro"
done < <(ls -1 *.fastq.gz | xargs -n 1 readlink | sed -r 's~/Samples_[^/]+/[^/]+$~~' | sort -u)
```

## Configuration

To resolve file(s) per sample a tab-separated **files.txt** file containing regular expression and sample basename can be added to the input or output directory:

### input/files.txt example
```
.*5610_89_CF25885_.*.gz 5610_89_CF25885
.*5610_8_CF13235_.*.gz  5610_8_CF13235
```

If not present the provided /input directory will be queried (recursively) for .fastq.gz files and output samplenames will consist of the the basename part before an index sequence 6 to 8bp in length.

To adjust pipeline parameters, use docker run with --env PARAM=value or --env-file=**config.txt** in the docker run command. See env_options.txt for parameters, their defaults and allowed alternative settings (variables therein are assigned as in a [bash] shell - don't use whitespace)

### config.txt example
```
MEM=3G
THREADS=4
PARALLEL=2

# breast or ovarian cancer
TYPE=breast # / ovarian

# b1.191 : breast BRCA1 191 probe classifier
# b1.371 : breast BRCA1 371 probe classifier
# b1 : for ovarian cancer BRCA1 classifier
# b2 : for breast or ovarian cancer BRCA2 classifier
BRCA_CLASS=b1.191 #/ b1.371 / b1 / b2

# CPFCOR: apply crossplatform correction
CPFCOR=TRUE # / FALSE

# SETM2C: set missing values to centroid mean (TRUE, null influence of the probe) or 
# linear interpolation of missing values (FALSE)
SETM2C=FALSE # / TRUE

# OUTCLS : output directory of classification. Note that trailing / is required for properly constructing
# the directory name: e.g. cls/

TAG=-GRCh38-blacklist-merged
```

======

## Running the container

### This runs the container  (full paths are required in docker)
docker run --rm -u $UID:$GROUPS \
  -v $path_to_output:/output:rw \
  $input \
  $reference \
  --env-file=config.txt \
  ovabrca:$tag

## ito execute manually or for development that does not require new software
docker run --rm -t -i -u $UID:$GROUPS \
  -v $path_to_output:/output:rw \
  $input \
  $reference \
  -v $path_to_app:/app:rw \
  --env-file=config.txt \
  ovabrca:$tag bash

## 
docker run --rm -u $UID:$GROUPS \
  -v $path/output/:/output:rw \
  $input \
  $reference \
  -v $another_path/GRCh38_51bp-q15-20k.bed.gz:/app/GRCh38_51bp-q15-20k.bed.gz:ro \
  --env PRODUCTION=FALSE \
  --env TAG=-test \
  --env SEQLEN=51 \
  ovabrca:0.9

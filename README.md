# Brcaness-classifier

To setup Docker use documentation [here](https://docs.docker.com/get-docker/) or provided by your Linux distribution.

## Building the docker image
Dependent on your docker configuration, you may have to add the flag '--network=host'.

```bash
tag=0.91
docker build --tag=ovabrca:$tag .
```

======
## Preparing docker container requirements

When the docker container is run, i.e. an instance of the built docker image, the [-v or --volume](https://docs.docker.com/engine/reference/run/#volume-shared-filesystems) option can be added to bind mount a directory (or a file) to a docker container. The docker container has no access to the local file system except for the configured volumes.

A docker container does abide to your system file permissions, so limit what you expose to a docker container. Read about security [here](https://docs.docker.com/engine/security/security/#docker-daemon-attack-surface).

While running, the brcaness-classifier docker requires /ref and /output volumes as well as /input, unless continuing with previous data.


### Prepare /ref volume

The /ref volume should contain a fasta that is bwa indexed; initialized for bwa. The docker image will try to index, if not already, but this require write access, time and a lot of RAM. One can choose to bwa index the fasta outside of the docker container, use bwa 0.7.17. To set up the alignment environment as used currently:

```bash
cd $genome_directory
wget ftp://ftp.ensembl.org/pub/release-76/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primaryassembly.fa.gz

## to bwa index outside of the container, run here:
# $path_to_bwa_0_7_17/bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa
# $path_to_samtools/samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

reference="-v `pwd`:/ref:rw"
cd -
```

## Input/Output access for the container

An input volume, containing raw data, as well as an empty output volume are required, when starting with a new set of samples.

```bash
mkdir -p output/ input/
input="-v $path_to_fastq_files:/input:ro"
```

The simple configuration is to place raw data in the provided $path_to_fastq_files directory. 

If intermediate files (cram files, counts or the classifier text files) are present, processing can be continued from there, providing the files pass consistency checks and are up-to-date. If no /input volume was provided, the files in the /output volume will be used as a starting point.

When using the input volume, mind that a docker container cannot follow symlinks. Files can be hard linked, if on the same physical disk. To include samples from multiple sources, multiple volumes can be nested under the /input volume. The $input variable can then be extended like this:

```bash
while read d; do
  p=input/$(basename $d)
  mkdir $p
  input="$input -v $d:/$p:ro"
done < <(ls -1 *.fastq.gz | xargs -n 1 readlink | sed -r 's~/Samples_[^/]+/[^/]+$~~' | sort -u)
```

The inclusion of nested input volumes are in particular required, if a same sample was sequenced in multiple runs, with an identical name. To combine the files in the pipeline the files.txt confiuration file (described further down below) requires regular expressions to match the source data to combine.

## Configuration

## Processing configuration

To adjust how samples are processed, environment variables in the docker are changed. Invoke docker run with --env PARAMETER=value or --env-file=**config.txt** in the docker run command. See env_options.txt for parameters that can be adjusted, their default values, and allowed alternative settings (variables therein are assigned as in a [bash] shell - don't use whitespace). Some parameters have a fourth (non-comment) column. Those can be set only if PRODUCTION=FALSE is set.

### config.txt example
```
MEM=3G
THREADS=4
PARALLEL=2

# breast or ovarian cancer
TYPE=breast

# b1.191 : breast BRCA1 191 probe classifier
# b1.371 : breast BRCA1 371 probe classifier
# b1 : for ovarian cancer BRCA1 classifier
# b2 : for breast or ovarian cancer BRCA2 classifier
BRCA_CLASS=b1.191

# CPFCOR: apply crossplatform correction
CPFCOR=TRUE

# SETM2C: set missing values to centroid mean (TRUE, null influence of the probe) or
# linear interpolation of missing values (FALSE)
SETM2C=FALSE

# OUTCLS : output directory of classification. Note that trailing / is required for properly constructing
# the directory name: e.g. cls/

TAG=-GRCh38-blacklist-merged
```

To configure multiple fastq file(s) for a single sample, a tab-separated **files.txt** file containing regular expression and sample output basename can be added to the /input or /output volume:

### input/files.txt example
```
.*5610_89_CF25885_.*.gz 5610_89_CF25885
.*5610_8_CF13235_.*.gz  5610_8_CF13235
```

If no files.txt is present in /input or /output volume, the provided /input volume will be queried, recursively, for .fastq.gz files and output samplenames will consist of the the basename part before an index sequence, 6 to 8bp in length. If no /input volume is present either, the /output volume will be quired for the classifier, NKI text file, count files, cram files to determine where to start.

======

## Running the container

This runs the container (full paths are required in docker). Make sure to set $tag to the correct docker image.
```
docker run --rm -u $UID:$GROUPS \
  -v $path_to_output:/output:rw \
  $input \
  $reference \
  --env-file=config.txt \
  ovabrca:$tag
```

To execute manually in a shell, or for development
```
docker run --rm -t -i -u $UID:$GROUPS \
  -v $path_to_output:/output:rw \
  $input \
  $reference \
  -v $path_to_app:/app:rw \
  --env-file=config.txt \
  ovabrca:$tag bash
```

51bp on GRCh38, workaround for old samples, requires mappability file.
```
docker run --rm -u $UID:$GROUPS \
  -v $path/output/:/output:rw \
  $input \
  $reference \
  -v $another_path/GRCh38_51bp-q15-20k.bed.gz:/app/GRCh38_51bp-q15-20k.bed.gz:ro \
  --env PRODUCTION=FALSE \
  --env TAG=-test \
  --env SEQLEN=51 \
  ovabrca:0.9
```

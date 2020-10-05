#!/usr/bin/env bash
# (c) Roel Kluin r.kluin@nki.nl - GPL v2

set -e              # any error is fatal..
set -o pipefail     # also if it occured in a pipe
set -o noclobber    # prevent file overwrites

warn() {
  echo -e "$1" 1>&2
}

die() {
  warn "$1"
  exit 1
}

# 1) check used environment variables, see env_option.txt
while read param default alt_regex affected comment; do
  if [ -z "${!param}" ]; then
    printf -v "$param" "$default"
  elif [ "${!param}" != "$default" ]; then
    warn "$param has value ${!param}, default is $default"
    [[ "${!param}" =~ ^$alt_regex$ ]] || die "Invalid value, should be $default or match /^$alt_regex$/"
    if [ -n "${affected}" -a "${affected:0:1}" != "#" ]; then
      warn "^^^^---------- modification of this parameter is EXPERIMENTAL, the $affected may be affected!"
      [ "$PRODUCTION" = TRUE ] && exit 1
    fi
  fi
done < <(sed -n -r '/^( *#.*)?$/b;p' /app/env_options.txt)

[ "$DRY_RUN" = "TRUE" ] && set -n # sort of dry run - tests for intermediate files do not resolve.
[ "$DEBUG" = "TRUE" ] && set -x   # print executed commands

# 2) make sure alignment prerequisites are present

# mount points available and correct permissions?
[ -d /input -a -r /input -a ! -w /input ]  || die "Required: -v \$path_to_input:/input:ro (check dir is readable)"
[ -d /output -a -r /output -a -w /output ] || die "Required: -v \$path_to_output:/output:rw (check dir is read-write)"

[ "$(ls -1a /output | wc -l)" -ne 2 ] && warn "Ouput dir is not empty"

# sufficient RAM?
got_mem=$(awk '{if ($1=="MemAvailable:") {print int($2/1048576)}}' /proc/meminfo)
test_available_mem() {
  (( $got_mem < ${1%G})) && die "Insufficient memory: ${got_mem}G (rounded down), I need $1"
}

test_available_mem $MEM || true # for actual alignment

# 3.1) check existence of output files, potentially skip to end
# final output files already present?
OUTCLS="/output/${OUTCLS}"

ratiosegm="/output/${CN}_1m${TAG}.txt"

do_classify() {
  mkdir -p $OUTCLS
  Rscript --verbose /app/Classification.R ${ratiosegm} ${TYPE} ${BRCA_CLASS} ${LEGACY} ${CPFCOR} ${SETM2C} ${OUTCLS}

  #Rscript /app/Classification.R \
  # --input ${ratiosegm} \
  # --tissue ${TYPE} \
  # --class ${BRCA_CLASS} \
  # --legacy FALSE \
  # --correct ${CPFCOR} \
  # --missing2centroid ${SETM2C} \
  # --outcls ${OUTCLS}
}

if [ -r $ratiosegm ]; then
  warn "Found $ratiosegm, skipping to classification. Remove if samples should be added, instead."
  do_classify
  exit 0
fi

# 3) otherwise prepare to run alignment, etc.


# 3.1) lock before any potential writing.
lock=/output/.brcaness_classifier.lock
trap 'rm -f "$lock"; exit $?' INT TERM EXIT KILL

touch $lock || die "no write access on /output directory"

bwaindex="/ref/$FASTA"

# Does the provided reference genome meet the pipeline requirements?
# Compare provided .fai with sequence lenghts expected (includes contig name check)
test_genome() {
  cmp <(cut -f -2 "${bwaindex}.fai") $CHROMINFO || die "Fasta differs in contig names, order or sequence lengths."

if [ -n $CHECKSUM ]; then
  # Test reference sequence. Comments or newline placements are ignored.
  echo "Verifying reference genome checksum.."
  sha256=$(sed '/^>/!b;s/[ \t].*$//;s/^/\t/' "${bwaindex}" | tr -d "\r\n" |
    sha256sum | cut -d " " -f 1)

  [ "$sha256" = $CHECKSUM ] || die "Mismatch in fasta sequence checksum"
fi
}

# 3.2) index genome if required
if [ ! -r "${bwaindex}" ]; then
  test_genome
  die "Required: -v \$path_to_ref:/ref:ro (containing gunzipped $FASTA)"
elif [ ! -r "${bwaindex}.bwt" -o ! -r "${bwaindex}.fai" ]; then
  [ -w /ref -a -r "${bwaindex}" ] || die "use -v \$path_to_ref:/ref:rw on first run"
  samtools faidx "${bwaindex}"
  test_genome
  warn "indexing reference genome - initialization - will take a long time.."
  test_available_mem 32G || true
  bwa index "${bwaindex}"
fi

# remove partial or outdated files.
rm_part() {
  for f in "$@"; do
    [ -f "$f" ] || continue;
    warn "Replacing outdated or broken file:\n$f"
    [ "$REMOVE_PARTIAL" = TRUE ] || die "Not replacing file(s) without --env REMOVE_PARTIAL=TRUE"
    [ "$f" -nt "$lock" ] && die "File $f was created twice in the pipeline. files.txt duplicates?"
    rm "$f"
  done
}


mkdir -p /output/log
[ -r "/input/files.txt" ] && files="/input/files.txt" || files="/output/files.txt"

bins=/app/windows-${KBIN_SIZE}k.txt
test_mapq_count() {
  [ $(wc -l < "$1") -ne $(wc -l < $bins) ] && return 2
  grep -w -q "[1-9][0-9]*$" "$1"
}

# loop over samples
if [ -r $files ]; then
  # a link does not work here
  cp $files /tmp/files.txt
  # pre check all files so we don't have to wait for an error.
  while read f etc; do
    [ -z "$f" ] && die "empty line in $files"
    files=($(find /input/ -type f -regex "^/input/${f#^}"))
    n=${#files[@]}
    [ $n -eq 0 ] && die "$files: no files for $f found"
    echo "found $n files for$f"

    for f in "${files[@]}"; do
      [ -z "$(zgrep -m 2 "^.*$" "$f" | sed -n -r "/^[ACTGN]{$SEQLEN}$/p")" ] &&
         die "$f is not sequence length $SEQLEN"
    done

    if [ "$f" != "${files[0]}" -a -z "$etc" ]; then
       die "$files: need basename column(2) for regex matched files."
    fi
  done < <(egrep -v "^(# .*)?$" /tmp/files.txt)

elif [ -d /input/ ]; then
  find /input/ -type f -name "*.fastq.gz" | sed 's~/input/~~' > /tmp/files.txt

  while read f; do
    [ -z "$(zgrep -m 2 "^.*$" "/input/$f" | sed -n -r "/^[ACTGN]{$SEQLEN}$/p")" ] &&
       die "$f is not sequence length $SEQLEN"
  done < <(cat /tmp/files.txt)
else
  find /output -type f -name "*.cram" | sed -r 's~/output/(.*)\.cram$~\1~' > /tmp/files.txt
fi

while read glob SM LB ID; do
 while [ $(jobs -p | wc -w) -ge $PARALLEL ]; do
  while read p; do
    # if signal cannot be sent process is finished, then get status code
    kill -0 $p 2> /dev/null || wait $p || die "Error ocurred: $?"
  done < <(jobs -p)
  sleep 1m
 done

 SM="${SM:-$(basename $glob | sed -r 's/_[ACTG]{6,}.*$//')}"

 LB="${LB:-$SM}"

 ID="${ID:-$(mktemp -u | cut -d "." -f2)}"

 RG="@RG\tID:$ID\tSM:$SM\tLB:$LB\tCN:$CN\tPL:$PL"

 mapqcount="/output/${LB}-counts-${KBIN_SIZE}000-q${MINQUAL}.txt"

 aln="/output/${LB}.cram"
 ai="${aln}.crai"

 if [ -f "$mapqcount" ]; then
   test_mapq_count "$mapqcount" && continue
   [ $? -eq 2 ] && warn "$mapqcount is truncated" || warn "$mapqcount: only zeroes"
   rm_part "$mapqcount"
 fi

 # per sample commands are chained and sent to background.

 # index is created after alignment. if missing alignment was incomplete
 if [ -e "$ai" ]; then
   # if quickcheck fails on existing cram file all steps are redone
   samtools quickcheck "$aln" || rm_part "$ai"
 fi &&

 if [ ! -f "$ai" ]; then

   rm_part "$aln" "/output/log/${LB}_alignment.log"
   warn "Aligning $glob to $aln"

   # samples can be split over multiple runs / requests.
   (find /input/ -type f -regex "^/input/${glob#^}" | xargs zcat |
   bwa mem -M -t $THREADS -R "$RG" "$bwaindex" -|
   samtools sort -m $MEM -@ 4 --reference "$bwaindex" \
       --output-fmt cram,version=3.0 \
       --output-fmt-option seqs_per_slice=100000 - -o "$aln" &&
   samtools index "$aln") &> "/output/log/${LB}_alignment.log"
   samtools quickcheck "$aln" || die "samtools quickcheck failed for $aln"

 fi &&

 warn "counting $aln into $mapqcount" &&

 samtools view --reference "$bwaindex" -bu -q $MINQUAL "$aln" |
 CRAM_REFERENCE="$bwaindex" bedtools coverage -counts -sorted \
   -g $CHROMINFO -a $bins -b stdin | bedtools sort > "$mapqcount" &&

 test_mapq_count "$mapqcount" || {
   [ $? -eq 2 ] && die "$mapqcount is truncated" || die "$mapqcount: only zeroes"
 }&

done < <(egrep -v "^(#.*)?$" /tmp/files.txt)

# wait for completion
wait 2> /dev/null || true

# create classifier file
cd /output
mkdir -p qc${KBIN_SIZE}K
[ -f "$ratiosegm" ] || Rscript /app/create_${CN}_1m.R ${KBIN_SIZE} ${MINQUAL} ${SEQLEN} ${BLACKLIST} $(basename "$ratiosegm" ".txt")

# run classifier
do_classify



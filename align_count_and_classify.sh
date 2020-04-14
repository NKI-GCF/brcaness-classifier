#!/bin/bash
# (c) Roel Kluin r.kluin@nki.nl - GPL v2

set -e
set -o noclobber


die() {
  echo -e "$1" 1>&2
  exit 1
}

# check that mount points are available and have right permissions
for d in input output ref; do
  [ -d /$d ] || die "Directory /$d is not mounted"
  [ -r "/$d" ] || die "/$d is not readable!"
  [ "/$d" != "/output" -a -w "$d" ] && die "$d is writable!"
done

# config mount is optional
if [ -d "/config" ]; then
  [ -r "/config" ] || die "config is not readable!"
  [ -w "/config" ] && die "config is writable!"
  if [ -r "/config/config.txt" ]; then
    for set_config in $(sed -n -r "s/^[ \t]*([A-Z_]{2,})[= \t]+([A-Za-z0-9._-]{1,})/\1=\2/p" /config/config.txt); do
      eval "$set_config"
    done
  fi

  # If a run was interupted and you need to continue where you left off
  # remove broken files, lockfile and set CONTINUE=true
  # samtools quickcheck and presence of .bai file may catch some errors
  [ "$(ls -1a /output | wc -l)" -ne 2 -a "$CONTINUE" != "true" ] && die "Ouput dir is not empty"

fi

# Defaults, do not set memory too high or this may crash.
[ -z "$PARALLEL" ] && PARALLEL=2
[ -z "$THREADS" ] && THREADS=4
[ -z "$MEM" ] && MEM=2G

# for readgroup. Not sure if the classifier is valid on a different platform
[ -z "$CN" ] && CN="NKI"
[ -z "$PL" ] && PL="Illumina"

# The filename can differ, checks below ensure the file is Ensembl GRCh37.75
[ -z "$FASTA" ] && FASTA="Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"

# These should not be changed for the classifier.
[ -z "$quality_cutoff" ] && quality_cutoff=37
[ -z "$kbin_size" ] && kbin_size=20

# check that all fasta and alignment prerequisites are present
bwaindex="/ref/$FASTA"
for f in "${bwaindex}"{,.bwt,.fai}; do

  [ -r "$f" ] || die "Cannot read file: $f"

done


lock=/output/.brcaness_classifier.lock
trap 'rm -f "$lock"; exit $?' INT TERM EXIT KILL

touch $lock || die "no write access on /output directory"

chromInfo=/app/ref/chromInfo.txt

ref_mismatch() {
  cp $chromInfo /output/
  die "Fasta .fai does not have the same order orlengths of contigs. Is it GrCh37.75?\n"\
      "See chromInfo.txt in output file for the expected contigs and lengths."
}

# Compare provided .fai with sequence lenghts expected
cmp <(cut -f -2 "${bwaindex}.fai") $chromInfo || ref_mismatch

if [ "$CHECKSUM" != "false" ]; then
  # Test reference sequence. Comments or newline placements are ignored.
  echo "Reference genome checksum.."
  sha256=$(cat "${bwaindex}" | sed '/^>/!b;s/[ \t].*$//;s/^/\t/' | tr -d "\r\n" |
    sha256sum | cut -d " " -f 1)

  [ "$sha256" = "fc1a03150e12c259cebe37aa6295ed612eee71cd643e529986697b7f13cc94c9" ] || ref_mismatch
fi

bins=/tmp/windows-${kbin_size}000.txt

bedtools makewindows -g "${bwaindex}.fai" -w ${kbin_size}000 > "${bins}"

mkdir /output/log


# loop over samples
if [ -r "/config/files.txt" ]; then
  cp /config/files.txt /tmp/files.txt
  # so we don't have to wait for the error.
  while read f etc; do
    [ -f "/input/$f" ] || die "files.txt:$f not found"
  done < <(egrep -v "^(# .*)?$" /tmp/files.txt)
else
  ls -1 /input/*/*.fastq.gz > /tmp/files.txt
fi
while read fastq ID SM LB; do

 if [ $(jobs | wc -l) -gt $PARALLEL ]; then
   sleep 60s
   continue
 fi

 [ -z "$ID" ] && ID="$(mktemp -u | cut -d "." -f2)"

 [ -z "$SM" ] && SM="$(basename "$fastq" | sed -r 's/_[ACTG]{6,}.*$//')"

 [ -z "$LB" ] && LB="$SM"

 mapqcount="${LB}-counts-${kbin_size}000-q${quality_cutoff}.txt"

 [ -e "$mapqcount" ] && continue

 echo -e "$fastq\t$mapqcount" >> /output/sample_list.txt

 bam="/output/${SM%.bam}.bam"

 {
  # align
  if [ ! -f "$bam.bai" -o -n "$(samtools quickcheck "$bam")" ]; then
   echo "Aligning $fastq to $bam" 1>&2
   [ -f "$mapqcount" ] && rm "$mapqcount"

   bwa mem -M -t $THREADS -R "@RG\tID:$ID\tSM:$SM\tLB:$LB\tCN:$CN\tPL:$PL" "$bwaindex" "$fastq" |
   samtools sort -m $MEM - -o "$bam" &&
   samtools index "$bam" &> "/output/log/${LB}_alignment.log"

  fi &&

  samtools quickcheck "$bam" &&

  if [ ! -f "$mapqcount" ]; then
   echo "counting $bam into $mapqcount" 1>&2 &&

   samtools view -q $quality_cutoff -bu "$bam" |
   bedtools coverage -counts -sorted -g $chromInfo -a "${bins}" -b stdin |
   bedtools sort > "$mapqcount" &> "/output/log/${LB}_mapq_counts.log"
   if [ $(cat "$mapqcount" | wc -l) -ne $(cat "$bins" | wc -l) ]; then
     echo "bincount nr of lines does not match mapqcounts" 1>&2
     kill $PID
   fi
  fi;
 }&
done < <(egrep -v "^(#.*)?$" /tmp/files.txt)

while [ $(jobs -p | wc -l) -ne 0 ]; do
 sleep 60s
done

# create classifier file
cd /output
Rscript /app/create_NKI_1m.R 20 2>&1 | tee /output/log/create_NKI_1m.log

# next: classifier
Rscript ClassifierStarter.R \
 -i "NKI_1M_output.txt"\
 -c "/app/database.json" \
 -l "/app" \
 -d "$DB_INSERT"


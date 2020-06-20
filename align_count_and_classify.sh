#!/bin/bash
# (c) Roel Kluin r.kluin@nki.nl - GPL v2

set -e
set -o noclobber
set -x


die() {
  echo -e "$1" 1>&2
  exit 1
}

# check that mount points are available and have correct permissions
for d in input ref; do
  [ -d /$d -a -r /$d -a ! -w /$d ] || die "Required: -v \$path_to_${d}:/$d:ro"
done

[ -d /output -a -r /output -a -w /output ] ||
    die "Required: -v \$path_to_output:/output:rw"

# config mount is optional
if [ -d /config ]; then
  [ -r /config -a ! -w /config ] ||
    die "Read-only required for config: -v \$path_to_config:/config:ro"

  # options are read from the config.txt file. eval!: only allow simple var assignments.
  if [ -r "/config/config.txt" ]; then
    for set_config in $(sed -n -r "s/^[ \t]*([A-Z][A-Za-z0-9_]+)[= \t]+([A-Za-z0-9._-]{1,})/\1=\2/p" /config/config.txt); do
      eval "$set_config"
    done
  fi

  # If a run was interupted and you need to continue where you left off
  # remove broken files, lockfile and set CONTINUE=true
  # the pipeline may also catch errors and replace corrupt files.
  [ "$(ls -1a /output | wc -l)" -ne 2 -a "$CONTINUE" != "true" ] && die "Ouput dir is not empty"

fi

[[ "$BRCA_NUM" =~ ^[12]$ ]] || die "BRCA_NUM missing in config.txt specify BRCA_NUM=1 or BRCA_NUM=2 (also for ovarium carcinoma)"

[[ "$TYPE" =~ ^[12]$ ]] || die "TYPE missing in config.txt TYPE=1 for mamma-, TYPE=2 for ovarium carcinoma."

# Defaults, do not set memory too high or this may crash.
[ -z "$PARALLEL" ] && PARALLEL=2
[ -z "$THREADS" ] && THREADS=4
[ -z "$MEM" ] && MEM=2G
[ -z "$SKIP_CRAM" ] && SKIP_CRAM=false

# for readgroup. Not sure if the classifier is valid on a different platform
[ -z "$CN" ] && CN="NKI"
[ -z "$PL" ] && PL="Illumina"

[ -z "$BUILD" ] && BUILD="GRCh38"

# The filename can differ, checks below ensure the file is Ensembl build
[ -z "$FASTA" ] && FASTA="Homo_sapiens.${BUILD}.dna.primary_assembly.fa"

# These should not be changed for the classifier.
[ -z "$QUALITY_CUTOFF" ] && QUALITY_CUTOFF=15
[ -z "$SEQUENCE_LENGTH" ] && SEQUENCE_LENGTH=65
[ -z "$KBIN_SIZE" ] && KBIN_SIZE=20



# check that all fasta and alignment prerequisites are present
bwaindex="/ref/$FASTA"
for f in "${bwaindex}"{,.bwt,.fai}; do

  [ -r "$f" ] || die "Cannot read file: $f"

done


lock=/output/.brcaness_classifier.lock
trap 'rm -f "$lock"; exit $?' INT TERM EXIT KILL

touch $lock || die "no write access on /output directory"

# noclobber prefents overwrites. With CONTINUE=true in config.txt overwrites are allowed, but
# only for files created before the lock file. If a file is created after then this was due to
# a duplicate in the files.txt list.
rm_old() {
  [ "$CONTINUE" = "true" ] || die "Overwriting file(s) only allowed with CONTINUE=true:\n$@"
  for f in "$@"; do
    [ -f "$f" ] || continue;
    [ "$f" -nt "$lock" ] && die "File $f was created twice in the pipeline (files.txt has duplicates?)"
    rm "$f"
  done
}

out="/output/NKI_1m${TAG}.txt"
outxlsx="${out%.txt}.xlsx"

chromInfo=/app/chromInfo.txt

ref_mismatch() {
  cp $chromInfo /output/
  die "Fasta .fai does not have the same order or lengths of contigs. Is it ${BUILD}?\n"\
      "See chromInfo.txt in output file for the expected contigs and lengths."
}

# Compare provided .fai with sequence lenghts expected
cmp <(cut -f -2 "${bwaindex}.fai") $chromInfo || ref_mismatch

if [ "$CHECKSUM" != "false" ]; then
  # Test reference sequence. Comments or newline placements are ignored.
  echo "Verifying reference genome checksum.."
  sha256=$(sed '/^>/!b;s/[ \t].*$//;s/^/\t/' "${bwaindex}" | tr -d "\r\n" |
    sha256sum | cut -d " " -f 1)

  [ "$sha256" = "257604babc88d1a24bffa13f41c39172681c000a1dd396b32ebdcde0d1fa78f6" ] || ref_mismatch
fi

bins=windows-${KBIN_SIZE}000.txt
[ -f "/input/$bins" ] && cp "/input/$bins" "/tmp/$bins" ||
  bedtools makewindows -g "${bwaindex}.fai" -w ${KBIN_SIZE}000 > "/tmp/${bins}"

gccontent="gccontent-${KBIN_SIZE}000.txt"
[ -f "/input/$gccontent" ] && cp "/input/$gccontent" "/tmp/$gccontent" ||
  bedtools nuc -fi "${bwaindex}" -bed "/tmp/${bins}" > "/tmp/$gccontent"

mkdir -p /output/log


# loop over samples
if [ -r "/config/files.txt" ]; then
  cp /config/files.txt /tmp/files.txt
  # so we don't have to wait for the error.
  while read f etc; do
    [ -z "$f" ] && die "empty line in files.txt"
    files=($(find /input/ -type f -regex "^/input/${f#^}"))

    n=${#files[@]}
    [ $n -eq 0 ] && die "files.txt: no files for $f not found"
    echo "found $n files for $f: ${files[@]}"

    if [ "/input/$f" != "${files[0]}" -a -z "$etc" ]; then
       die "files.txt: need basename column(3) for regex matched files."
    fi
  done < <(egrep -v "^(# .*)?$" /tmp/files.txt)
else
  find /input/ -type f -name "*.fastq.gz" | sed 's~/input/~~' > /tmp/files.txt
fi

now_running=0
(
while read f SM LB ID; do
 now_running=$(((now_running+1)%(PARALLEL+1)))
 [ $now_running -eq 0 ] && wait

 [ -z "$SM" ] && SM="$(basename $f | sed -r 's/_[ACTG]{6,}.*$//')"

 [ -z "$LB" ] && LB="$SM"

 [ -z "$ID" ] && ID="$(mktemp -u | cut -d "." -f2)"

 mapqcount="/output/${LB}-counts-${KBIN_SIZE}000-q${QUALITY_CUTOFF}.txt"
 aln="/output/${LB}.cram"
 [ ! -e "${aln}.crai" ] && aln="/output/${LB}.bam"
 ai="$(echo ${aln} | sed -n -r 's/(^.+\.(b|cr)*am)$/\1.\2ai/p')"

 {
  # index is created after alignment.
  while [ -e "$ai" ]; do
    samtools quickcheck "$aln" && break;
    rm_old "$aln" "$ai"
    if [[ $aln =~ \.cram$ ]]; then
      aln="${aln%.cram}.bam"
      ai="${ai%.cram.crai}.bam.bai"
      [ -e "$ai" ] || rm_old "$mapqcount" "/output/log/${LB}_mapq_counts.log" "$out"
    else
     # if quickcheck fails on existing bam file all steps are redone
      rm_old "$mapqcount" "/output/log/${LB}_mapq_counts.log" "$out"
    fi
  done

  if [ ! -f "$ai" ]; then
   echo "Aligning $f to $aln" 1>&2
   # noclobber is active: ensure mapqcount can be done
   [ -f "$mapqcount" ] && rm_old "$mapqcount" "/output/log/${LB}_mapq_counts.log"

   # samples can be split over multiple runs / requests.
   (find /input/ -type f -regex "^/input/${f#^}" | xargs zcat |
   bwa mem -M -t $THREADS -R "@RG\tID:$ID\tSM:$SM\tLB:$LB\tCN:$CN\tPL:$PL" "$bwaindex" -|
   samtools sort -m $MEM -@ 4 - -o "$aln" &&
   samtools index "$aln") &> "/output/log/${LB}_alignment.log"
   samtools quickcheck "$aln" || die "samtools quickcheck failed for $aln"
  fi

  if [ ! -f "$mapqcount" ]; then
   echo "counting $aln into $mapqcount" 1>&2 &&

   samtools view --reference "$bwaindex" -bu -q $QUALITY_CUTOFF "$aln" |
   CRAM_REFERENCE="$bwaindex" bedtools coverage -counts -sorted \
    -g $chromInfo -a "/tmp/${bins}" -b stdin | bedtools sort > "$mapqcount"

   if [ $(wc -l < "$mapqcount") -ne $(wc -l < "/tmp/$bins") ]; then
     die "$bincount nr of lines does not match mapqcounts"
   else
     grep -w -q "[1-9][0-9]*$" "$mapqcount" || die "$bincount contains only zeros"
   fi
  fi
  if [ "$SKIP_CRAM" != "true" ]; then
    if [[ $aln =~ \.bam$ ]]; then
      # conversion to cram after counts: order for the classifier
      cram="${aln%.bam}.cram"
      samtools view --reference "$bwaindex" \
       --output-fmt cram,version=3.0 \
       --output-fmt-option seqs_per_slice=100000 $aln -o "$cram"
      samtools index "$cram"
      samtools quickcheck "$cram" || die "samtools quickcheck [bam -> cram] failed for $aln"
      rm "$aln" "$ai"
    fi
  fi
 }&
done < <(egrep -v "^(#.*)?$" /tmp/files.txt)

# wait for completion
for i in seq 0 $PARALLEL; do
  wait 2> /dev/null || true
done
)

# create classifier file
cd /output
mkdir -p qc${KBIN_SIZE}K
[ -f "$out" ] || Rscript /app/create_NKI_1m.R ${KBIN_SIZE} ${QUALITY_CUTOFF} ${SEQUENCE_LENGTH} ${BLACKLIST} NKI_1m${TAG} 2>&1 | tee /output/log/create_NKI_1m${TAG}.log

#[ -f "$outxlsx" ] || Rscript /app/cmd_EL.R "$out"

# next: classifier
#Rscript /app/classifierR.R \
# -i "$outxlsx" \
# -s $(wc -l < /tmp/files.txt) \
# -t "$TYPE" \
# -o "/output/NKI_1M_${TYPE}_classifier" \
# -b "$BRCA_NUM"


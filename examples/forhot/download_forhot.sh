#!/usr/bin/env bash
# Download the 20 ForHot runs (SRA PRJNA663238) from ENA over HTTPS, check
# every file against ENA's MD5, then join the runs of each sample into
# fastq/<sample>_1/_2.fastq.gz (the names the samplesheet expects).
# About 232 GB of compressed reads, plus joined copies for the three samples
# sequenced in several runs (single-run samples are symlinked).
#
# Usage: bash examples/forhot/download_forhot.sh /path/to/forhot
#        (then set fastq_dir to /path/to/forhot/fastq in the config)
set -euo pipefail
export LC_ALL=C   # md5sum must print "OK"/"FAILED", whatever the system language
HERE=$(cd "$(dirname "$0")" && pwd)
TARGET=${1:-.}
mkdir -p "$TARGET" && cd "$TARGET"
mkdir -p raw fastq metadata
RUNS=metadata/runs_forhot.tsv
tr -d '\r' < "$HERE/runs_forhot.tsv" > "$RUNS"   # tolerate Windows line endings

JOBS=${JOBS:-4}
echo "== download over HTTPS ($JOBS files at a time; wget -c resumes partial files, skips complete ones)"
if ! tail -n +2 "$RUNS" | awk -F'\t' '{print $6; print $8}' | sed 's#^ftp://#https://#' \
     | xargs -n 1 -P "$JOBS" wget -c -q --tries=20 --waitretry=60 -P raw; then
  echo "Some downloads failed: run the script again, it resumes where it stopped." >&2
  exit 1
fi

echo "== md5 check"
tail -n +2 "$RUNS" | awk -F'\t' '{print $7"  raw/"$3"_1.fastq.gz"; print $9"  raw/"$3"_2.fastq.gz"}' > metadata/ena_md5.txt
CHK=metadata/md5_check_$(date +%F).txt
md5sum -c metadata/ena_md5.txt > "$CHK" 2>&1 || true
cat "$CHK"
n_exp=$(wc -l < metadata/ena_md5.txt); n_ok=$(grep -c ': OK$' "$CHK" || true)
if [ "$n_ok" -ne "$n_exp" ]; then echo "MD5 CHECK FAILED: $n_ok of $n_exp files OK. Stop." >&2; exit 1; fi
echo "MD5: all $n_ok files OK"
chmod a-w raw/*.fastq.gz

echo "== join runs per sample (same run order for both mates)"
for s in $(tail -n +2 "$RUNS" | cut -f1 | sort -u); do
  runs=$(awk -F'\t' -v s="$s" '$1==s' "$RUNS" | sort -t$'\t' -k2,2n | cut -f3)
  for m in 1 2; do
    out=fastq/${s}_${m}.fastq.gz
    [ -s "$out" ] && continue
    files=$(for r in $runs; do printf 'raw/%s_%s.fastq.gz ' "$r" "$m"; done)
    n=$(echo $runs | wc -w)
    if [ "$n" -eq 1 ]; then ln -s "$(readlink -f $files)" "$out"
    else cat $files > "$out"; fi
    echo "$out <- $files"
  done
done
date > metadata/download_done.txt
echo "done"

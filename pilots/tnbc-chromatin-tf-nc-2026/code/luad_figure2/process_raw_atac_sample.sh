#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "Usage: process_raw_atac_sample.sh <run_root> <system> <sample_id> <sample_slug> <run> <layout>" >&2
  exit 2
fi

RUN_ROOT=$(realpath "$1")
SYSTEM=$2
SAMPLE_ID=$3
SLUG=$4
RUN=$5
LAYOUT=$6

case "$SYSTEM" in PDX|cell_line) ;; *) echo "Invalid system: $SYSTEM" >&2; exit 2;; esac
case "$LAYOUT" in PAIRED|SINGLE) ;; *) echo "Invalid layout: $LAYOUT" >&2; exit 2;; esac
[[ "$SLUG" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "Unsafe sample slug: $SLUG" >&2; exit 2; }

REF="$RUN_ROOT/reference"
RAW="$RUN_ROOT/data/raw/sra"
OUT="$RUN_ROOT/data/processed/raw_atac/$SYSTEM/$SLUG"
TMP="$RUN_ROOT/tmp_processing/$SYSTEM/$SLUG"
LOG="$RUN_ROOT/logs/raw_atac/$SYSTEM"
QC="$RUN_ROOT/audit/raw_atac_qc/$SYSTEM"
CUTADAPT="$RUN_ROOT/tools/cutadapt-venv/bin/cutadapt"
MACS2="$RUN_ROOT/tools/macs2-venv/bin/macs2"
BLACKLIST="$REF/hg38-blacklist.v2.bed"
mkdir -p "$RAW" "$OUT" "$TMP" "$LOG" "$QC"

for required in "$CUTADAPT" "$MACS2" "$BLACKLIST"; do
  if [[ ! -e "$required" ]]; then echo "Missing required reference/tool: $required" >&2; exit 2; fi
done
if [[ ! -s "$REF/hg38_bowtie2/hg38.1.bt2" && ! -s "$REF/hg38_bowtie2/hg38.1.bt2l" ]]; then
  echo "Missing human Bowtie2 index" >&2
  exit 2
fi
if [[ "$SYSTEM" == "PDX" && ! -s "$REF/mm10_bowtie2/mm10.1.bt2" && ! -s "$REF/mm10_bowtie2/mm10.1.bt2l" ]]; then
  echo "Missing mouse Bowtie2 index" >&2
  exit 2
fi

COMPLETE="$OUT/.complete"
if [[ -s "$COMPLETE" ]]; then
  echo "SKIP_COMPLETE $SYSTEM $SAMPLE_ID $RUN"
  exit 0
fi

exec > >(tee "$LOG/${SLUG}.log") 2>&1
echo "START_UTC=$(date -u +%FT%TZ) system=$SYSTEM sample=$SAMPLE_ID run=$RUN layout=$LAYOUT"

prefetch --max-size 100G --output-directory "$RAW" "$RUN"
SRA_PATH="$RAW/$RUN/$RUN.sra"
if [[ ! -s "$SRA_PATH" ]]; then
  SRA_PATH=$(find "$RAW/$RUN" -maxdepth 2 -type f -name '*.sra' -print -quit)
fi
[[ -s "$SRA_PATH" ]] || { echo "SRA file not found for $RUN" >&2; exit 1; }

mkdir -p "$TMP/fasterq_tmp" "$TMP/fastq"
if [[ "$LAYOUT" == "PAIRED" ]]; then
  if [[ ! -s "$TMP/fastq/${RUN}_1.fastq" || ! -s "$TMP/fastq/${RUN}_2.fastq" ]]; then
    rm -f "$TMP/fastq/${RUN}_1.fastq" "$TMP/fastq/${RUN}_2.fastq"
    fasterq-dump --threads 4 --split-files --temp "$TMP/fasterq_tmp" --outdir "$TMP/fastq" "$SRA_PATH"
  fi
else
  if [[ ! -s "$TMP/fastq/${RUN}.fastq" && ! -s "$TMP/fastq/${RUN}_1.fastq" ]]; then
    fasterq-dump --threads 4 --split-files --temp "$TMP/fasterq_tmp" --outdir "$TMP/fastq" "$SRA_PATH"
  fi
fi

ADAPTER="CTGTCTCTTATACACATCT"
if [[ "$LAYOUT" == "PAIRED" ]]; then
  R1="$TMP/fastq/${RUN}_1.fastq"
  R2="$TMP/fastq/${RUN}_2.fastq"
  [[ -s "$R1" && -s "$R2" ]] || { echo "Missing paired FASTQ output" >&2; exit 1; }
  "$CUTADAPT" -j 4 -a "$ADAPTER" -A "$ADAPTER" -m 20:20 --pair-filter=any \
    -o "$TMP/trimmed_R1.fastq.gz" -p "$TMP/trimmed_R2.fastq.gz" "$R1" "$R2" \
    > "$LOG/${SLUG}.cutadapt.log"
  bowtie2 --very-sensitive --no-mixed --no-discordant -p 8 \
    -x "$REF/mm10_bowtie2/mm10" \
    -1 "$TMP/trimmed_R1.fastq.gz" -2 "$TMP/trimmed_R2.fastq.gz" \
    --un-conc-gz "$TMP/nonmouse_%.fastq.gz" -S /dev/null \
    2> "$LOG/${SLUG}.mouse_depletion.log"
  bowtie2 --very-sensitive --no-mixed --no-discordant -p 8 \
    -x "$REF/hg38_bowtie2/hg38" \
    -1 "$TMP/nonmouse_1.fastq.gz" -2 "$TMP/nonmouse_2.fastq.gz" \
    2> "$LOG/${SLUG}.human_alignment.log" \
    | samtools view -b -o "$TMP/aligned.bam" -
else
  R1="$TMP/fastq/${RUN}.fastq"
  if [[ ! -s "$R1" ]]; then R1="$TMP/fastq/${RUN}_1.fastq"; fi
  [[ -s "$R1" ]] || { echo "Missing single-end FASTQ output" >&2; exit 1; }
  "$CUTADAPT" -j 4 -a "$ADAPTER" -m 20 \
    -o "$TMP/trimmed_R1.fastq.gz" "$R1" > "$LOG/${SLUG}.cutadapt.log"
  bowtie2 --very-sensitive -p 8 -x "$REF/hg38_bowtie2/hg38" \
    -U "$TMP/trimmed_R1.fastq.gz" 2> "$LOG/${SLUG}.human_alignment.log" \
    | samtools view -b -o "$TMP/aligned.bam" -
fi

samtools sort -n -@ 4 -o "$TMP/name_sorted.bam" "$TMP/aligned.bam"
samtools fixmate -m -@ 4 "$TMP/name_sorted.bam" "$TMP/fixmate.bam"
samtools sort -@ 4 -o "$TMP/coordinate_sorted.bam" "$TMP/fixmate.bam"
samtools markdup -r -s -@ 4 "$TMP/coordinate_sorted.bam" "$TMP/deduplicated.bam" \
  2> "$LOG/${SLUG}.markdup.log"

if [[ "$LAYOUT" == "PAIRED" ]]; then
  samtools view -b -f 2 -F 2820 -q 30 "$TMP/deduplicated.bam" \
    | samtools view -h - \
    | awk 'BEGIN{OFS="\t"} /^@/{print; next} $3!="chrM" && $3!="chrY"{print}' \
    | samtools view -b -o "$TMP/qc_preblacklist.bam" -
  # pairToBed requires query-name grouped mates.  Keep or discard the entire
  # fragment so a blacklisted mate cannot leave an apparent singleton.
  samtools sort -n -@ 4 -o "$TMP/qc_preblacklist_name_sorted.bam" "$TMP/qc_preblacklist.bam"
  bedtools pairtobed -abam "$TMP/qc_preblacklist_name_sorted.bam" -b "$BLACKLIST" -type neither \
    > "$TMP/blacklist_filtered.bam"
else
  samtools view -b -F 2820 -q 30 "$TMP/deduplicated.bam" \
    | samtools view -h - \
    | awk 'BEGIN{OFS="\t"} /^@/{print; next} $3!="chrM" && $3!="chrY"{print}' \
    | samtools view -b -o "$TMP/qc_preblacklist.bam" -
  bedtools intersect -v -abam "$TMP/qc_preblacklist.bam" -b "$BLACKLIST" \
    > "$TMP/blacklist_filtered.bam"
fi
samtools sort -@ 4 -o "$OUT/${SLUG}.filtered.bam" "$TMP/blacklist_filtered.bam"
samtools index "$OUT/${SLUG}.filtered.bam"

mkdir -p "$OUT/peaks"
if [[ "$LAYOUT" == "PAIRED" ]]; then
  "$MACS2" callpeak -t "$OUT/${SLUG}.filtered.bam" -f BAMPE -g hs \
    --keep-dup all -B --nomodel --SPMR -q 0.01 \
    -n "$SLUG" --outdir "$OUT/peaks" > "$LOG/${SLUG}.macs2.log" 2>&1
  FILTERED_READS=$(samtools view -c "$OUT/${SLUG}.filtered.bam")
  FILTERED_UNITS=$((FILTERED_READS / 2))
else
  bedtools bamtobed -i "$OUT/${SLUG}.filtered.bam" \
    | awk 'BEGIN{OFS="\t"} {if($6=="+"){$2+=4;$3+=4}else{$2-=5;$3-=5;if($2<0){$2=0}} print}' \
    > "$OUT/${SLUG}.tn5_shifted.bed"
  "$MACS2" callpeak -t "$OUT/${SLUG}.tn5_shifted.bed" -f BED -g hs \
    --keep-dup all -B --shift -75 --extsize 150 --nomodel --SPMR -q 0.01 \
    -n "$SLUG" --outdir "$OUT/peaks" > "$LOG/${SLUG}.macs2.log" 2>&1
  FILTERED_READS=$(samtools view -c "$OUT/${SLUG}.filtered.bam")
  FILTERED_UNITS=$FILTERED_READS
fi

PEAKS="$OUT/peaks/${SLUG}_peaks.narrowPeak"
[[ -s "$PEAKS" ]] || { echo "MACS2 produced no narrowPeak for $SAMPLE_ID" >&2; exit 1; }
PEAK_COUNT=$(wc -l < "$PEAKS")
IN_PEAK_READS=$(bedtools intersect -u -abam "$OUT/${SLUG}.filtered.bam" -b "$PEAKS" | samtools view -c -)
FRIP=$(awk -v a="$IN_PEAK_READS" -v b="$FILTERED_READS" 'BEGIN{if(b>0) printf "%.8f",a/b; else print "NA"}')
QC_STATUS=PASS
if (( FILTERED_UNITS < 1000000 || PEAK_COUNT < 10000 )); then QC_STATUS=FAIL; fi

printf 'system\tsample_id\tsample_slug\trun\tlibrary_layout\tfiltered_human_reads\tfiltered_human_units\tpeak_count\tfrip\thard_qc_status\n' \
  > "$QC/${SLUG}.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$SYSTEM" "$SAMPLE_ID" "$SLUG" "$RUN" "$LAYOUT" "$FILTERED_READS" "$FILTERED_UNITS" "$PEAK_COUNT" "$FRIP" "$QC_STATUS" \
  >> "$QC/${SLUG}.tsv"

printf 'completed_utc=%s\nqc_status=%s\n' "$(date -u +%FT%TZ)" "$QC_STATUS" > "$COMPLETE"
echo "COMPLETE system=$SYSTEM sample=$SAMPLE_ID units=$FILTERED_UNITS peaks=$PEAK_COUNT frip=$FRIP qc=$QC_STATUS"

# Raw .sra, final filtered BAM, shifted BED, peaks and all logs are retained.
# Only mechanically regenerable per-sample scratch files are removed.
rm -f \
  "$TMP/aligned.bam" "$TMP/name_sorted.bam" "$TMP/fixmate.bam" \
  "$TMP/coordinate_sorted.bam" "$TMP/deduplicated.bam" \
  "$TMP/qc_preblacklist.bam" "$TMP/qc_preblacklist_name_sorted.bam" \
  "$TMP/blacklist_filtered.bam"

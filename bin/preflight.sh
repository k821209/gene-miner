#!/usr/bin/env bash
# preflight.sh -- check, before any task is submitted, the things that otherwise
# fail hours into a run.
#
# Why this exists. Every expensive failure in this pipeline so far has had the
# same shape: something cheap to verify was wrong, nothing looked wrong, and the
# process that needed it runs late. Four, measured, in one week:
#
#   * conda env prefixes pointing at $HOME/miniconda3 on a host whose envs live
#     elsewhere -> BUSCO_TRAIN dies on `cp: cannot stat .../augustus/config`,
#     after every cached stage has been re-checked.
#   * the eggNOG database absent -> EGGNOG fails after BUILD_UNION, i.e. after
#     masking, alignment, assembly and BOTH ab-initio streams are paid for.
#   * eggnog-mapper installed against python 3.12+, where distutils is gone ->
#     emapper.py cannot start at all. Same late failure, and the install looks
#     clean.
#   * a host unable to reach busco-data.ezlab.org -> BUSCO retries for ~17 min
#     and then fails the run, 17 h in.
#
# None of these need a genome, a read, or a core to detect. So they are checked
# here, in seconds, and the run stops with a message that names the fix.
#
# Three verdicts, not two. FAIL stops the run. WARN is a real risk that is not
# certainly fatal. UNKNOWN is "this check could not run", which is reported
# separately and loudly, because a check that cannot say it was skipped is worse
# than no check: it reads as a pass.
#
# Usage: preflight.sh --genome F --reads GLOB --proteome F --busco-lineage X \
#          --env-annot D --env-augustus D --env-rmod D --env-busco D \
#          --env-eggnog D --env-genemark D --genemark-etp-dir D --eggnog-db D \
#          --run-genemark 0|1 --run-qc 0|1 [--augustus-species S] \
#          [--masked-genome F] [--masked-rmout F] [--repeat-lib F] [--repeat-gff F]
#          [--ref-gff F]
set -uo pipefail

GENOME= READS= PROTEOME= LINEAGE= AUG_SP=
E_ANNOT= E_AUG= E_RMOD= E_BUSCO= E_EGG= E_GM= GM_ETP= EGGDB=
RUN_GM=0 RUN_QC=1 MASKED_G= MASKED_O= REPEAT_LIB= REPEAT_GFF= REF_GFF=
while [ $# -gt 0 ]; do
  case $1 in
    --genome) GENOME=$2; shift 2;;
    --reads) READS=$2; shift 2;;
    --proteome) PROTEOME=$2; shift 2;;
    --busco-lineage) LINEAGE=$2; shift 2;;
    --augustus-species) AUG_SP=$2; shift 2;;
    --env-annot) E_ANNOT=$2; shift 2;;
    --env-augustus) E_AUG=$2; shift 2;;
    --env-rmod) E_RMOD=$2; shift 2;;
    --env-busco) E_BUSCO=$2; shift 2;;
    --env-eggnog) E_EGG=$2; shift 2;;
    --env-genemark) E_GM=$2; shift 2;;
    --genemark-etp-dir) GM_ETP=$2; shift 2;;
    --eggnog-db) EGGDB=$2; shift 2;;
    --run-genemark) RUN_GM=$2; shift 2;;
    --run-qc) RUN_QC=$2; shift 2;;
    --masked-genome) MASKED_G=$2; shift 2;;
    --masked-rmout) MASKED_O=$2; shift 2;;
    --repeat-lib) REPEAT_LIB=$2; shift 2;;
    --repeat-gff) REPEAT_GFF=$2; shift 2;;
    --ref-gff) REF_GFF=$2; shift 2;;
    *) echo "preflight: unknown option $1" >&2; exit 2;;
  esac
done

NFAIL=0 NWARN=0 NUNK=0
pass(){ printf '  ok      %s\n' "$*"; }
fail(){ printf '  FAIL    %s\n' "$*"; NFAIL=$((NFAIL+1)); }
warn(){ printf '  WARN    %s\n' "$*"; NWARN=$((NWARN+1)); }
unk(){  printf '  UNKNOWN %s\n' "$*"; NUNK=$((NUNK+1)); }
skip(){ printf '  skip    %s\n' "$*"; }

istrue(){ case "${1,,}" in 1|true|yes|y|on) return 0;; *) return 1;; esac; }

# ---------------------------------------------------------------- inputs
echo "[preflight] inputs"
if [ -z "$GENOME" ]; then fail "--genome not given"
elif [ ! -s "$GENOME" ]; then fail "genome not found or empty: $GENOME"
else
  if [ "$(head -c1 "$GENOME")" != ">" ]; then
    fail "genome does not start with '>': $GENOME (gzipped? pass an uncompressed FASTA)"
  else
    NSEQ=$(grep -c '^>' "$GENOME")
    pass "genome: $NSEQ seqid(s), $(du -h "$GENOME" | cut -f1)"
    # RepeatMasker could not handle this project's 50-character lettuce seqids
    # (Dovetail_09Sept_Map_inspected_12-07-2015_1_v8_lg_1) and the run failed
    # AFTER RepeatModeler had spent 9 h. The exact threshold was not pinned
    # down, so this warns with the measured case rather than asserting a limit.
    MAXID=$(grep '^>' "$GENOME" | sed 's/^>//' | awk '{print length($1)}' | sort -rn | head -1)
    LONGEST=$(grep '^>' "$GENOME" | sed 's/^>//' | awk '{print length($1), $1}' | sort -rn | head -1 | cut -d' ' -f2-)
    if [ "${MAXID:-0}" -ge 50 ] && [ -z "$MASKED_G" ]; then
      warn "longest seqid is $MAXID chars ($LONGEST).
             RepeatMasker rejected 50-character names in this project and the run
             failed after RepeatModeler had already run for 9 h. If MASK_GENOME
             fails on sequence names, shorten them (bin/shorten_seqids.py if
             present, keeping the map) and restore them before comparing against
             a reference annotation."
    else
      pass "longest seqid: ${MAXID:-?} chars"
    fi
    # An unwrapped FASTA (one sequence per line) fails deep inside BUSCO: its
    # bbtools step runs at a fixed --xmx=120m and dies with a Java
    # OutOfMemoryError that names neither the file nor the reason. Measured here
    # on a genome whose longest line was 37,915,939 characters; rewrapping to 60
    # columns fixed it. Cheap to see, impossible to diagnose from the symptom.
    MAXLINE=$(awk '!/^>/ {if (length($0) > m) m = length($0)} END {print m+0}' "$GENOME")
    if [ "$MAXLINE" -gt 1000000 ]; then
      warn "longest sequence line is $MAXLINE characters -- this genome is
             effectively unwrapped. BUSCO's bbtools step runs at a fixed
             --xmx=120m and dies with a Java OutOfMemoryError naming nothing.
             Rewrap first, e.g.  seqkit seq -w 60  (or any 60-80 column rewrap)."
    else
      pass "sequence lines wrapped (longest $MAXLINE chars)"
    fi
    if grep -qm1 '^>lcl|' "$GENOME"; then
      warn "genome headers carry NCBI 'lcl|' prefixes. If the annotation GFF
             does not, nothing will join and every gene looks novel -- strip
             them first (this cost a full re-run here)."
    fi
  fi
fi

if [ -z "$PROTEOME" ]; then fail "--proteome not given"
elif [ ! -s "$PROTEOME" ]; then fail "proteome not found or empty: $PROTEOME"
elif [ "$(head -c1 "$PROTEOME")" != ">" ]; then fail "proteome is not FASTA: $PROTEOME"
else pass "proteome: $(grep -c '^>' "$PROTEOME") sequence(s)"; fi

if [ -z "$READS" ]; then
  fail "--reads not given"
else
  # The glob is the pipeline's own {1,2} form; count mate-1 files via shell globbing
  # of the pattern with the brace replaced, which is what Channel.fromFilePairs does.
  R1=$(eval ls -1 "${READS//\{1,2\}/1}" 2>/dev/null | wc -l)
  R2=$(eval ls -1 "${READS//\{1,2\}/2}" 2>/dev/null | wc -l)
  if [ "$R1" -eq 0 ]; then
    fail "--reads matched no files: $READS"
  elif [ "$R1" -ne "$R2" ]; then
    fail "read pairs are unbalanced: $R1 mate-1 file(s) vs $R2 mate-2 (glob: $READS)"
  else
    pass "reads: $R1 pair(s)"
  fi
fi

# ---------------------------------------------------------------- masked injection
echo "[preflight] masking"
if [ -n "$REPEAT_GFF" ] && [ -n "$MASKED_G" ]; then
  fail "--repeat_gff and --masked_genome are both set; they are two ways to do the
             same thing. Pass one."
fi
if [ -n "$REPEAT_GFF" ] && [ -z "$MASKED_G" ]; then
  if [ ! -s "$REPEAT_GFF" ]; then fail "--repeat_gff not found or empty: $REPEAT_GFF"
  else
    # Checked here because the conversion's one success-shaped failure is a GFF
    # whose seqids do not match the genome: it yields an empty .out, which a
    # repeat-poor genome also yields legitimately.
    GFFSEQ=$(awk -F'\t' '!/^#/ && NF>=9 {print $1}' "$REPEAT_GFF" | sort -u)
    if [ -n "$GENOME" ] && [ -s "$GENOME" ]; then
      FASEQ=$(grep '^>' "$GENOME" | sed 's/^>//' | awk '{print $1}' | sort -u)
      SHARED=$(comm -12 <(printf '%s\n' "$GFFSEQ") <(printf '%s\n' "$FASEQ") | grep -c . || true)
      NGFF=$(printf '%s\n' "$GFFSEQ" | grep -c . || true)
      if [ "${SHARED:-0}" -eq 0 ]; then
        fail "--repeat_gff shares NO seqid with the genome ($NGFF distinct in the GFF).
             The conversion would produce an empty .out, which looks like a
             repeat-poor genome rather than a mismatch. Check for prefixes
             ('lcl|', 'chr') or a different assembly."
      else
        pass "repeat GFF: $(awk -F'\t' '!/^#/ && NF>=9' "$REPEAT_GFF" | wc -l) interval(s), $SHARED/$NGFF seqid(s) shared with the genome"
      fi
    else
      unk "cannot compare --repeat_gff seqids against the genome (genome unreadable)"
    fi
  fi
  skip "RepeatMasker/RepeatModeler (masking derived from --repeat_gff, then gated by VALIDATE_MASKED_INPUT)"
elif [ -n "$MASKED_G" ]; then
  [ -s "$MASKED_G" ] || fail "--masked_genome not found: $MASKED_G"
  if [ -z "$MASKED_O" ]; then fail "--masked_genome given without --masked_rmout (TE_FILTER needs the coordinates)"
  elif [ ! -s "$MASKED_O" ]; then fail "--masked_rmout not found: $MASKED_O"
  else pass "injected masking: $(basename "$MASKED_G") + $(basename "$MASKED_O") (VALIDATE_MASKED_INPUT will compare them)"; fi
  skip "RepeatMasker/RepeatModeler (masking is injected)"
else
  if [ -z "$E_RMOD" ] || [ ! -d "$E_RMOD" ]; then
    fail "masking will run but the rmod env is missing: ${E_RMOD:-<unset>}
             set GM_CONDA_BASE to the prefix holding envs/, or --env_rmod directly."
  else
    if PATH="$E_RMOD/bin:$PATH" command -v RepeatMasker >/dev/null 2>&1; then
      pass "RepeatMasker found in $E_RMOD"
      # RepeatMasker >= 4.1.5 aborts at startup without a FamDB, even with
      # --repeat_lib, so an empty famdb is fatal rather than a reduced run.
      FAMDIR=$(ls -d "$E_RMOD"/share/RepeatMasker/Libraries/famdb 2>/dev/null || true)
      if [ -z "$FAMDIR" ]; then
        fail "no Libraries/famdb under $E_RMOD/share/RepeatMasker -- RepeatMasker
             will not start (it aborts even when --repeat_lib is supplied)."
      else
        NH5=$(ls "$FAMDIR"/*.h5 2>/dev/null | wc -l)
        if [ "$NH5" -eq 0 ]; then
          fail "famdb at $FAMDIR holds no .h5 partitions -- RepeatMasker reports
             Families:0 and dies. Fetch them (setup_envs.sh does)."
        else
          pass "famdb: $NH5 partition(s)"
        fi
      fi
    else
      fail "RepeatMasker not executable in $E_RMOD/bin"
    fi
  fi
  if [ -n "$REPEAT_LIB" ]; then
    if [ -s "$REPEAT_LIB" ]; then pass "repeat library: $(grep -c '^>' "$REPEAT_LIB") families (de-novo RepeatModeler skipped)"
    else fail "--repeat_lib not found or empty: $REPEAT_LIB"; fi
  fi
fi

# ---------------------------------------------------------------- conda envs
echo "[preflight] conda envs"
# Checked by RUNNING a binary, not by looking for a directory: an env can exist
# and still be unusable, which is how the eggnog env got through.
check_env(){ # name prefix cmd...
  local label=$1 prefix=$2; shift 2
  if [ -z "$prefix" ]; then fail "$label env prefix is unset"; return; fi
  if [ ! -d "$prefix" ]; then
    fail "$label env not found: $prefix
             nextflow.config resolves prefixes from GM_CONDA_BASE (default
             \$HOME/miniconda3). Export GM_CONDA_BASE=<prefix holding envs/>."
    return
  fi
  local missing=""
  for c in "$@"; do PATH="$prefix/bin:$PATH" command -v "$c" >/dev/null 2>&1 || missing="$missing $c"; done
  if [ -n "$missing" ]; then fail "$label env at $prefix is missing:$missing"
  else pass "$label env: $(basename "$prefix") ($*)"; fi
}
check_env annot    "$E_ANNOT"  hisat2 stringtie samtools gffread
# TransDecoder is NOT in the env's bin/: the TRANSDECODER process puts
# opt/transdecoder/util on PATH and takes the three entry points from there.
# Checking bin/ instead reported a healthy install as broken, which would have
# blocked every run -- so this mirrors what the process actually does.
if [ -n "$E_ANNOT" ] && [ -d "$E_ANNOT" ]; then
  TDU="$E_ANNOT/opt/transdecoder/util"
  TDMISS=""
  for t in TransDecoder.LongOrfs TransDecoder.Predict cdna_alignment_orf_to_genome_orf.pl; do
    [ -x "$TDU/$t" ] || TDMISS="$TDMISS $t"
  done
  if [ -n "$TDMISS" ]; then fail "TransDecoder entry points missing from $TDU:$TDMISS"
  else pass "transdecoder util: 3/3 in opt/transdecoder/util"; fi
fi
check_env augustus "$E_AUG"    augustus
if [ -n "$E_AUG" ] && [ ! -d "$E_AUG/config" ]; then
  fail "no config/ in the augustus env ($E_AUG) -- BUSCO_TRAIN copies it and
             fails with 'cp: cannot stat'. This is the symptom of a wrong
             GM_CONDA_BASE."
fi
check_env busco    "$E_BUSCO"  busco

# ---------------------------------------------------------------- BUSCO lineage
echo "[preflight] BUSCO lineage"
if [ -z "$LINEAGE" ]; then
  fail "--busco_lineage not given"
elif [ -d "$LINEAGE" ]; then
  if [ "$(basename "$(dirname "$LINEAGE")")" = "lineages" ]; then
    NHMM=$(ls "$LINEAGE"/hmms 2>/dev/null | wc -l)
    if [ "$NHMM" -eq 0 ]; then fail "staged lineage $LINEAGE has no hmms/ -- incomplete download"
    else
      pass "staged lineage: $(basename "$LINEAGE"), $NHMM HMM(s) -- will run --offline, no network needed"
      CFGNAME=$(sed -n 's/^name=//p' "$LINEAGE/dataset.cfg" 2>/dev/null | head -1)
      if [ -n "$CFGNAME" ] && [ "$CFGNAME" != "$(basename "$LINEAGE")" ]; then
        warn "the dataset's own dataset.cfg says name=$CFGNAME while the directory
             is $(basename "$LINEAGE"). BUSCO's distribution does this (solanales_odb10
             ships name=solanaceae_odb10); the directory name is what --offline
             uses. Worth a sentence in any report that prints BUSCO's output."
      fi
    fi
  else
    fail "$LINEAGE is a directory but not inside a lineages/ dir; --offline needs
             <download_path>/lineages/<name>."
  fi
else
  # A bare name means BUSCO will fetch it, and it contacts ezlab on EVERY start
  # even when the dataset is already cached. One host here could not reach ezlab
  # while two siblings in the same zone could.
  if command -v curl >/dev/null 2>&1; then
    if curl -sfL --max-time 25 -o /dev/null https://busco-data.ezlab.org/v5/data/file_versions.tsv 2>/dev/null; then
      pass "lineage '$LINEAGE' will be downloaded; busco-data.ezlab.org is reachable"
    else
      warn "lineage '$LINEAGE' will be downloaded, but busco-data.ezlab.org is NOT
             reachable from this host. BUSCO retries ~17 min and then fails the run.
             Stage the dataset and pass the DIRECTORY instead
             (<dir>/lineages/$LINEAGE), which switches BUSCO to --offline."
    fi
  else
    unk "cannot test ezlab reachability (no curl); lineage '$LINEAGE' will be downloaded"
  fi
fi

# ---------------------------------------------------------------- GeneMark-ETP
echo "[preflight] GeneMark-ETP (3rd stream)"
if ! istrue "$RUN_GM"; then
  skip "GENEMARK_ETP (--run_genemark is off)"
else
  check_env genemark "$E_GM" perl python3
  if [ -z "$GM_ETP" ] || [ ! -x "$GM_ETP/bin/gmetp.pl" ]; then
    fail "no gmetp.pl at ${GM_ETP:-<unset>}/bin/gmetp.pl -- clone GeneMark-ETP there
             or set GENEMARK_ETP_DIR."
  else
    pass "gmetp.pl: $GM_ETP/bin/gmetp.pl"
    MISSING_T=""
    for t in bedtools samtools hisat2 diamond stringtie gffread; do
      [ -x "$GM_ETP/tools/$t" ] || MISSING_T="$MISSING_T $t"
    done
    if [ -n "$MISSING_T" ]; then fail "GeneMark-ETP tools/ is missing:$MISSING_T (run_genemark_etp.sh puts tools/ first on PATH)"
    else pass "GeneMark-ETP bundled tools: 6/6"; fi
    if [ -n "$E_GM" ] && [ -x "$E_GM/bin/perl" ]; then
      MISSING_M=""
      for m in YAML YAML::XS Math::Utils Statistics::LineFit Parallel::ForkManager Hash::Merge MCE List::MoreUtils; do
        PATH="$E_GM/bin:$PATH" "$E_GM/bin/perl" -M"$m" -e1 >/dev/null 2>&1 || MISSING_M="$MISSING_M $m"
      done
      if [ -n "$MISSING_M" ]; then fail "genemark env perl is missing:$MISSING_M (cpanm them into the env)"
      else pass "genemark perl modules: 8/8"; fi
    else
      unk "cannot test perl modules (no perl in ${E_GM:-<unset>}/bin)"
    fi
  fi
fi

# ---------------------------------------------------------------- QC stage
echo "[preflight] QC stage (eggNOG / taxonomy / TE / BUSCO)"
if ! istrue "$RUN_QC"; then
  skip "EGGNOG and the QC filters (--run_qc is off)"
else
  check_env eggnog "$E_EGG" emapper.py diamond
  # emapper.py must RUN, not merely exist: installed against python 3.12+ it
  # cannot import distutils and dies on startup, while the install looks clean
  # and the failure waits until after BUILD_UNION.
  if [ -n "$E_EGG" ] && [ -x "$E_EGG/bin/emapper.py" ]; then
    if EMOUT=$(PATH="$E_EGG/bin:$PATH" "$E_EGG/bin/emapper.py" --version 2>&1); then
      pass "emapper runs: $(printf '%s' "$EMOUT" | grep -oE 'emapper-[0-9.]+' | head -1)"
    else
      if printf '%s' "$EMOUT" | grep -q "No module named 'distutils'"; then
        fail "emapper.py cannot start: no distutils. The eggnog env has python
             $(PATH="$E_EGG/bin:$PATH" "$E_EGG/bin/python" -c 'import sys;print("%d.%d"%sys.version_info[:2])' 2>/dev/null || echo '?'),
             and python >= 3.12 removed it while eggnog-mapper still imports it.
             Fix: $E_EGG/bin/python -m pip install setuptools
             (or recreate the env with \"python<3.12\")."
      else
        fail "emapper.py --version failed: $(printf '%s' "$EMOUT" | tail -1)"
      fi
    fi
  fi
  if [ -z "$EGGDB" ]; then
    fail "--eggnog-db is unset; EGGNOG runs after BUILD_UNION and will fail there."
  else
    DBMISS=""
    for f in eggnog.db eggnog_proteins.dmnd; do [ -s "$EGGDB/$f" ] || DBMISS="$DBMISS $f"; done
    if [ -n "$DBMISS" ]; then
      fail "eggNOG DB at $EGGDB is missing:$DBMISS
             EGGNOG runs AFTER masking, alignment, assembly, both ab-initio
             streams and BUILD_UNION, so this fails at the most expensive moment.
             Fetch it before starting (bin/run_eggnog.sh documents the archives)."
    else
      HDR=$(head -c 15 "$EGGDB/eggnog.db" 2>/dev/null)
      if [ "$HDR" != "SQLite format 3" ]; then
        fail "$EGGDB/eggnog.db is not a SQLite file (starts '$HDR') -- a truncated
             or half-written download has the right name and the wrong contents."
      else
        pass "eggNOG DB: eggnog.db $(du -h "$EGGDB/eggnog.db" | cut -f1) (SQLite ok), dmnd $(du -h "$EGGDB/eggnog_proteins.dmnd" | cut -f1)"
      fi
    fi
  fi
fi

# ---------------------------------------------------------------- reference to supplement
echo "[preflight] reference annotation to supplement (--ref_gff)"
if [ -z "$REF_GFF" ]; then
  skip "MERGE_REFERENCE (--ref_gff not given; only our own catalogue is produced)"
elif [ ! -s "$REF_GFF" ]; then
  fail "--ref_gff not found or empty: $REF_GFF"
else
  case "$REF_GFF" in
    *.gz) fail "--ref_gff is gzipped: $REF_GFF
             gm_compare and the merge read plain text. Decompress it first." ;;
    *)
      NG=$(awk -F'\t' '!/^#/ && $3=="gene"' "$REF_GFF" | wc -l)
      if [ "$NG" -eq 0 ]; then
        fail "--ref_gff has no 'gene' rows in column 3: $REF_GFF
             the merge keys on gene records; an annotation using another feature
             type would silently merge nothing."
      else
        pass "reference annotation: $NG gene record(s)"
        # The failure this prevents: mismatched sequence names share nothing, so
        # every locus of ours is classed novel and the merge appends the whole
        # catalogue instead of supplementing. No error is raised either way.
        if [ -n "$GENOME" ] && [ -s "$GENOME" ]; then
          RS=$(awk -F'\t' '!/^#/ && NF>=9 {print $1}' "$REF_GFF" | sort -u)
          GS=$(grep '^>' "$GENOME" | sed 's/^>//' | awk '{print $1}' | sort -u)
          SH=$(comm -12 <(printf '%s\n' "$RS") <(printf '%s\n' "$GS") | grep -c . || true)
          NR=$(printf '%s\n' "$RS" | grep -c . || true)
          if [ "${SH:-0}" -eq 0 ]; then
            fail "--ref_gff shares NO seqid with the genome ($NR distinct in the GFF).
             Every locus would be reported novel and the merge would append the
             whole catalogue rather than supplement it. Check for prefixes
             ('lcl|', 'chr') or shortened names needing restoration first."
          else
            pass "reference seqids shared with the genome: $SH/$NR"
          fi
        else
          unk "cannot compare --ref_gff seqids against the genome (genome unreadable)"
        fi
      fi ;;
  esac
fi

# ---------------------------------------------------------------- verdict
echo
if [ "$NFAIL" -gt 0 ]; then
  echo "[preflight] $NFAIL problem(s) that will stop this run, $NWARN warning(s), $NUNK unchecked."
  echo "[preflight] Nothing has been submitted. Fix the FAIL lines and re-run."
  exit 1
fi
if [ "$NUNK" -gt 0 ]; then
  echo "[preflight] $NUNK check(s) could NOT run (UNKNOWN above) -- that is not a pass."
fi
echo "[preflight] ok${NWARN:+ ($NWARN warning(s))}${NUNK:+, $NUNK unchecked}"
exit 0

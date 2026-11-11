#!/usr/bin/env bash
# write_provenance.sh -- record WHICH versions produced this output.
#
# "Annotated with Gene-Miner" is not a reproducible statement. The gene count
# and the TE filtering move with the RepeatMasker build and the Dfam release
# (4.1.5 + Dfam 3.7 and 4.2.4 + Dfam 4.0 are different experiments, and the
# second has 7.6 M families against the first's 19,768), with the AUGUSTUS model
# source, and with the eggNOG DB version. A methods section needs those named,
# and a reader asking "which Dfam?" six months later has no other way to find
# out.
#
# Every lookup is defensive and an unknown is PRINTED as unknown rather than
# omitted: a version that silently disappears from this file reads as "not
# applicable", which is the one thing it must not read as.
#
# Usage: write_provenance.sh <out_file> <final_gff> [<busco_summary_dir>]
set -uo pipefail
OUT=${1:?usage: write_provenance.sh out_file final_gff [busco_dir]}
GFF=${2:-}
BUSCODIR=${3:-}

UNKNOWN=0
# 2>&1, not 2>/dev/null: several of these tools print their version to stderr
# (augustus does), and discarding it turned a determinable version into
# "unknown".
v(){ # v <label> <regex-or-empty> <command...>
  local label=$1 pat=$2; shift 2
  local val
  val=$("$@" 2>&1 | tr -d '\r')
  if [ -n "$pat" ]; then val=$(printf '%s' "$val" | grep -oE "$pat" | head -1); else val=$(printf '%s' "$val" | head -1); fi
  lit "$label" "$val"
}
# An EMPTY value must not be printed as an empty field. It happened: a parse that
# returned nothing printed "Dfam families:" followed by whitespace, which reads
# as a value of none rather than as a failed lookup. So emptiness is funnelled
# through one place.
lit(){
  if [ -z "${2:-}" ]; then printf '  %-26s unknown (lookup returned nothing)\n' "$1:"; UNKNOWN=$((UNKNOWN+1));
  else printf '  %-26s %s\n' "$1:" "$2"; fi
}
unk(){ printf '  %-26s unknown (%s)\n' "$1:" "$2"; UNKNOWN=$((UNKNOWN+1)); }
# Version from the conda metadata, for a tool with no --version flag
# (TransDecoder has none -- it exits with "don't understand options").
conda_ver(){ # conda_ver <env-prefix> <package-name>
  local f
  f=$(ls "$1"/conda-meta/"$2"-*.json 2>/dev/null | head -1) || true
  [ -n "$f" ] && basename "$f" .json | sed "s/^$2-//" || true
}

E_ANNOT=${GM_PV_ANNOT:-}; E_AUG=${GM_PV_AUGUSTUS:-}; E_RMOD=${GM_PV_RMOD:-}
E_BUSCO=${GM_PV_BUSCO:-}; E_EGG=${GM_PV_EGGNOG:-}; GM_ETP=${GM_PV_GENEMARK_ETP:-}
EGGDB=${GM_PV_EGGNOG_DB:-}; LINEAGE=${GM_PV_LINEAGE:-}; PROJDIR=${GM_PV_PROJECTDIR:-}

{
echo "Gene-Miner provenance"
echo "generated: $(date -u +%Y-%m-%dT%H:%M:%SZ) on $(hostname)"
echo
echo "pipeline"
if [ -n "$PROJDIR" ] && [ -d "$PROJDIR/.git" ]; then
  lit "gene-miner commit" "$(git -C "$PROJDIR" log -1 --format='%H %ci' 2>/dev/null || echo unknown)"
  lit "gene-miner branch" "$(git -C "$PROJDIR" rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
  DIRTY=$(git -C "$PROJDIR" status --porcelain 2>/dev/null | wc -l)
  [ "${DIRTY:-0}" -eq 0 ] && lit "working tree" "clean" || lit "working tree" "MODIFIED ($DIRTY file(s)) -- the commit above does not fully describe this run"
else
  unk "gene-miner commit" "no git checkout at ${PROJDIR:-<unset>}"
fi
v "nextflow" 'version [0-9.]+( build [0-9]+)?' nextflow -version

echo
echo "repeat masking  (changes the gene count and the TE filter)"
if [ -n "$E_RMOD" ] && [ -d "$E_RMOD" ]; then
  export PATH="$E_RMOD/bin:$PATH"
  v "RepeatMasker" '[0-9]+\.[0-9.]+' RepeatMasker -v
  v "RepeatModeler" '[0-9]+\.[0-9.]+' bash -c 'RepeatModeler -h 2>&1 | grep -m1 -oE "RepeatModeler - [0-9.]+"'
  v "rmblastn" '[0-9]+\.[0-9.]+\+?' rmblastn -version
  FAMDB=$(ls -d "$E_RMOD"/share/RepeatMasker/Libraries/famdb 2>/dev/null)
  if [ -n "$FAMDB" ]; then
    FPY=$(ls "$E_RMOD"/share/famdb-*/famdb.py 2>/dev/null | head -1)
    if [ -n "$FPY" ]; then
      INFO=$(python3 "$FPY" -i "$FAMDB" info 2>&1)
      lit "Dfam version" "$(printf '%s' "$INFO" | sed -n 's/^ *Version *: *//p' | head -1)"
      lit "Dfam release date" "$(printf '%s' "$INFO" | sed -n '/^Database/,$p' | sed -n 's/^ *Date *: *//p' | head -1)"
      lit "famdb format" "$(printf '%s' "$INFO" | sed -n 's/^ *FamDB Creation Format Version *: *//p' | head -1)"
      lit "famdb partitions" "$(ls "$FAMDB"/*.h5 2>/dev/null | wc -l)"
      # famdb 3.x `info` reports partitions, not a family total; the count comes
      # from a per-clade query, and asking for one here would mean guessing the
      # clade. Named as unavailable rather than left blank.
      unk "Dfam family count" "famdb $(printf '%s' "$INFO" | sed -n 's/^ *FamDB Creation Format Version *: *//p' | head -1) info reports partitions only -- use famdb.py info <clade>"
    else
      unk "Dfam version" "famdb.py not found under $E_RMOD/share"
    fi
  else
    unk "Dfam version" "no Libraries/famdb (masking may have been injected)"
  fi
else
  lit "repeat masking" "not run from this env (injected --masked_genome/--repeat_gff, or env unset)"
fi

echo
echo "gene prediction"
if [ -n "$E_AUG" ] && [ -d "$E_AUG" ]; then PATH="$E_AUG/bin:$PATH" v "AUGUSTUS" '[0-9]+\.[0-9.]+' augustus --version
else unk "AUGUSTUS" "env unset"; fi
lit "AUGUSTUS model" "${GM_PV_AUG_MODEL:-unknown}"
if [ -n "$GM_ETP" ] && [ -d "$GM_ETP/.git" ]; then
  lit "GeneMark-ETP commit" "$(git -C "$GM_ETP" log -1 --format='%h %ci' 2>/dev/null || echo unknown)"
elif [ -n "$GM_ETP" ] && [ -d "$GM_ETP" ]; then
  lit "GeneMark-ETP" "$GM_ETP (no git metadata)"
else
  lit "GeneMark-ETP" "not used (--run_genemark off)"
fi

echo
echo "transcript evidence"
if [ -n "$E_ANNOT" ] && [ -d "$E_ANNOT" ]; then
  export PATH="$E_ANNOT/bin:$PATH"
  v "HISAT2" 'version [0-9.]+' hisat2 --version
  v "StringTie" '[0-9]+\.[0-9.]+' stringtie --version
  v "samtools" '[0-9]+\.[0-9.]+' bash -c "samtools --version"
  TDU="$E_ANNOT/opt/transdecoder/util"
  if [ -x "$TDU/TransDecoder.LongOrfs" ]; then
    TDV=$(conda_ver "$E_ANNOT" transdecoder)
    if [ -n "$TDV" ]; then lit "TransDecoder" "$TDV (conda metadata; the scripts have no --version)"
    else unk "TransDecoder" "present at $TDU but no conda metadata to read a version from"; fi
  else unk "TransDecoder" "not found in $TDU"; fi
else unk "transcript tools" "annot env unset"; fi

echo
echo "functional annotation / QC"
if [ -n "$E_EGG" ] && [ -d "$E_EGG" ]; then
  EMV=$(PATH="$E_EGG/bin:$PATH" "$E_EGG/bin/emapper.py" --version 2>&1)
  lit "eggnog-mapper" "$(printf '%s' "$EMV" | grep -oE 'emapper-[0-9.]+' | head -1 || echo unknown)"
  lit "eggNOG DB expected" "$(printf '%s' "$EMV" | grep -oE 'Expected eggNOG DB version: [0-9.]+' | sed 's/.*: //' | head -1 || echo unknown)"
  PATH="$E_EGG/bin:$PATH" v "diamond" '[0-9]+\.[0-9.]+' bash -c "diamond --version"
else lit "eggnog-mapper" "not used (--run_qc off)"; fi
if [ -n "$EGGDB" ] && [ -s "$EGGDB/eggnog.db" ]; then
  lit "eggNOG DB file" "$EGGDB/eggnog.db ($(stat -c%s "$EGGDB/eggnog.db") bytes, mtime $(stat -c%y "$EGGDB/eggnog.db" | cut -d' ' -f1))"
else lit "eggNOG DB file" "not used"; fi
if [ -n "$E_BUSCO" ] && [ -d "$E_BUSCO" ]; then PATH="$E_BUSCO/bin:$PATH" v "BUSCO" '[0-9]+\.[0-9.]+' busco --version
else unk "BUSCO" "env unset"; fi
# The lineage DIRECTORY name is what BUSCO uses; the dataset's own dataset.cfg can
# disagree with it (solanales_odb10 ships name=solanaceae_odb10), so record both.
if [ -n "$LINEAGE" ]; then
  lit "BUSCO lineage (given)" "$LINEAGE"
  if [ -d "$LINEAGE" ] && [ -r "$LINEAGE/dataset.cfg" ]; then
    lit "  dataset.cfg name" "$(sed -n 's/^name=//p' "$LINEAGE/dataset.cfg" | head -1)"
    lit "  dataset.cfg created" "$(sed -n 's/^creation_date=//p' "$LINEAGE/dataset.cfg" | head -1)"
    lit "  BUSCOs in dataset" "$(sed -n 's/^number_of_BUSCOs=//p' "$LINEAGE/dataset.cfg" | head -1)"
  fi
else unk "BUSCO lineage" "unset"; fi
if [ -n "$BUSCODIR" ] && [ -d "$BUSCODIR" ]; then
  S=$(ls "$BUSCODIR"/short_summary*.txt 2>/dev/null | head -1)
  [ -n "$S" ] && lit "BUSCO result" "$(grep -oE 'C:[0-9.]+%\[S:[0-9.]+%,D:[0-9.]+%\],F:[0-9.]+%,M:[0-9.]+%,n:[0-9]+' "$S" | head -1 || echo 'summary present, line not matched')"
fi

echo
echo "inputs"
for pair in "genome:${GM_PV_GENOME:-}" "proteome:${GM_PV_PROTEOME:-}" "repeat_lib:${GM_PV_REPEAT_LIB:-}" \
            "masked_genome:${GM_PV_MASKED_GENOME:-}" "repeat_gff:${GM_PV_REPEAT_GFF:-}"; do
  k=${pair%%:*}; f=${pair#*:}
  if [ -z "$f" ]; then continue; fi
  if [ -s "$f" ]; then lit "$k" "$f ($(stat -c%s "$f") bytes, mtime $(stat -c%y "$f" | cut -d' ' -f1))"
  else unk "$k" "$f not readable now"; fi
done
lit "reads" "${GM_PV_READS:-unknown}"

echo
echo "output"
if [ -n "$GFF" ] && [ -s "$GFF" ]; then
  lit "final gff" "$GFF"
  lit "genes" "$(awk -F'\t' '$3=="gene"' "$GFF" | wc -l)"
  lit "mRNAs" "$(awk -F'\t' '$3=="mRNA"' "$GFF" | wc -l)"
else unk "final gff" "${GFF:-<unset>} not readable"; fi

echo
if [ "$UNKNOWN" -gt 0 ]; then
  echo "NOTE: $UNKNOWN field(s) could not be determined and are marked 'unknown' above."
  echo "      They are listed rather than omitted: an absent line would read as"
  echo "      'not applicable', which is not what it means."
else
  echo "All fields determined."
fi
} | tee "$OUT"
exit 0

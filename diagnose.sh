#!/usr/bin/env bash
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT" || exit 1

echo "================ 1. CODE MODULES ================"
for f in \
  plastanno/__init__.py plastanno/pipeline.py \
  plastanno/core/feature.py plastanno/core/reconcile.py \
  plastanno/identify/ir_detector.py plastanno/identify/closest_rel.py \
  plastanno/identify/engine_a.py plastanno/identify/engine_b.py \
  plastanno/annotate/special_cases.py plastanno/output/writers.py \
  plastanno.py scripts/build/build_all.py scripts/build/build_protein_db.py ; do
  if [ -f "$f" ]; then printf "  %-44s %5s dòng\n" "$f" "$(wc -l < "$f")"
  else printf "  %-44s  THIẾU\n" "$f"; fi
done

echo ""
echo "================ 2. DATABASE ================"
for d in blast_db exon_db hmm_db protein_db rrna_db trna_db ; do
  if [ -d "database/$d" ]; then
    printf "  %-12s %6s  (%s tệp)\n" "$d" "$(du -sh database/$d 2>/dev/null|cut -f1)" "$(find database/$d -type f|wc -l)"
  else printf "  %-12s  THIẾU\n" "$d"; fi
done
echo "  protein_db (phải có .faa/.fasta): $(ls database/protein_db/ 2>/dev/null | head -3 | tr '\n' ' ')"
echo "  gene_catalog.json: $([ -f database/gene_catalog.json ] && echo CO || echo THIEU)"

echo ""
echo "================ 3. IMPORT TEST ================"
python3 - << 'PY'
import importlib
for m in ["plastanno.pipeline","plastanno.core.feature","plastanno.core.reconcile",
          "plastanno.identify.ir_detector","plastanno.identify.closest_rel",
          "plastanno.identify.engine_a","plastanno.identify.engine_b",
          "plastanno.annotate.special_cases","plastanno.output.writers"]:
    try: importlib.import_module(m); print("  OK  ", m)
    except Exception as e: print("  FAIL", m, "->", type(e).__name__, e)
PY

echo ""
echo "================ 4. EXTERNAL TOOLS ================"
for t in makeblastdb blastn tblastn exonerate nhmmer hmmsearch aragorn tRNAscan-SE barrnap mafft muscle ; do
  printf "  %-14s %s\n" "$t" "$(command -v $t 2>/dev/null || echo KHONG-TIM-THAY)"
done

echo ""
echo "================ 5. 5 TỆP .py SỬA GẦN NHẤT ================"
find . -name "*.py" -not -path "*__pycache__*" -printf "%TY-%Tm-%Td %TH:%TM  %p\n" 2>/dev/null | sort -r | head -5

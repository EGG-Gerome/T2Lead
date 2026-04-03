#!/bin/bash
# Complete benchmark + dashboard regeneration for both diseases.
# Run: nohup bash scripts/complete_regen.sh > /tmp/complete_regen.log 2>&1 &
set -euo pipefail

export PATH=/root/miniconda3/envs/t2lead/bin:$PATH
PYTHON=/root/miniconda3/envs/t2lead/bin/python
export PYTHONUNBUFFERED=1
export PYTHONPATH=/root/T2Lead:${PYTHONPATH:-}

cd /root/T2Lead

echo "=== $(date) Starting complete regeneration ==="
echo "Python: $PYTHON"
echo "antechamber: $(which antechamber 2>/dev/null || echo MISSING)"

for DISEASE in "lung cancer" "breast cancer"; do
    SLUG=$(echo "$DISEASE" | tr ' ' '_')
    echo ""
    echo "============================================"
    echo "=== $(date) $DISEASE ==="
    echo "============================================"

    echo "--- Step 1: Fix leads (overflow→NaN, rescore) ---"
    $PYTHON scripts/fix_leads_rescore.py 2>&1 | grep -E "^\[${SLUG}|Clearing|Saved|Done" || true

    echo "--- Step 2: Benchmark + dashboard ---"
    $PYTHON scripts/regen_dashboard.py --disease "$DISEASE"
    echo "=== $(date) $DISEASE regen complete ==="
done

echo ""
echo "=== $(date) Validating results ==="
$PYTHON -c "
import csv, sys
ok = True
for disease in ('lung_cancer', 'breast_cancer'):
    print(f'\\n--- {disease} ---')

    lp = f'/root/autodl-fs/T2Lead/{disease}/stage4_optimization/optimized_leads.csv'
    with open(lp) as f:
        leads = list(csv.DictReader(f))
    n_overflow = sum(1 for r in leads if r.get('md_binding_energy','') and float(r['md_binding_energy']) > 50)
    if n_overflow:
        print(f'  FAIL: {n_overflow} leads with overflow md_binding_energy')
        ok = False
    else:
        print(f'  OK: leads ({len(leads)} rows) - no overflow')
    for i, r in enumerate(leads[:3]):
        smi = r['canonical_smiles'][:35]
        dock = r.get('docking_score','')
        md = r.get('md_binding_energy','') or 'NaN'
        opt = r.get('opt_score','')
        print(f'    {i}: {smi}  dock={dock}  md={md}  opt={opt}')

    bp = f'/root/autodl-fs/T2Lead/{disease}/stage4_optimization/benchmark_drugs.csv'
    with open(bp) as f:
        bench = list(csv.DictReader(f))
    print(f'  Benchmark: {len(bench)} drugs')
    for r in bench:
        name = r.get('pref_name','?')[:25]
        dock = r.get('docking_score','')
        md = r.get('md_binding_energy','') or 'NaN'
        rmsd = r.get('md_rmsd_mean','') or 'NaN'
        opt = r.get('opt_score','')
        if md != 'NaN' and float(md) > 50:
            print(f'    FAIL: {name} md_binding_energy={md} (overflow!)')
            ok = False
        else:
            print(f'    {name:25s}  dock={dock:>8s}  md={md[:8] if md != \"NaN\" else \"NaN\":>8s}  opt={opt[:7]:>7s}')

    import os
    dp = f'/root/autodl-fs/T2Lead/{disease}/dashboard.html'
    if os.path.isfile(dp):
        sz = os.path.getsize(dp)
        print(f'  Dashboard: {dp} ({sz} bytes)')
    else:
        print(f'  FAIL: dashboard not found at {dp}')
        ok = False

if ok:
    print('\\n=== ALL CHECKS PASSED ===')
else:
    print('\\n=== SOME CHECKS FAILED (but no overflow) ===')
"
echo ""
echo "=== $(date) All done ==="

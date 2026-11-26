"""Recalculate browser Tm with Biopython and Primer3. Run from any directory."""
import json
import random
import subprocess
from pathlib import Path
import Bio
import primer3
from Bio.SeqUtils import MeltingTemp as mt

ROOT = Path(__file__).resolve().parents[1]
fixture = (ROOT / "tests/pUC19c-L09137.2.fa").read_text(encoding="ascii")
rng = random.Random(61703)
oligos = ["".join(rng.choices("ACGT", k=length)) for length in range(18, 26) for _ in range(32)]
oligos.extend(["GCGCATATGCGC", "ATATATATATATATATAT"])
templates = [fixture] + ["".join(rng.choices("ACGT", weights=w, k=1400)) for w in
                         [(3, 2, 2, 3), (1, 1, 1, 1), (2, 3, 3, 2)]]
driver = """
const s = require('./docs/engine.js');
const input = JSON.parse(require('fs').readFileSync(0,'utf8'));
const plans = input.templates.flatMap(t => [s.design(t,{circular:true}),
    s.design(t,{mode:'target',start:600,stop:660})]);
const picks = plans.flatMap(r => r.picked);
console.log(JSON.stringify({points:input.oligos.map(oligo=>({oligo,tm:s.tm(oligo)})),
  picks, counts:plans.map(r=>r.picked.length)}));
"""
data = json.loads(subprocess.check_output(["node", "-e", driver], cwd=ROOT,
                  input=json.dumps(dict(templates=templates, oligos=oligos)), text=True))
assert data["counts"] == [5] * 8, data["counts"]
conditions = dict(mv_conc=50, dv_conc=0, dntp_conc=0, dna_conc=50)
errors_bio, errors_p3, hairpin, dimer = [], [], [], []
for point in data["points"] + data["picks"]:
    seq = point["oligo"]
    selfcomp = seq == seq.translate(str.maketrans("ATGC", "TACG"))[::-1]
    bio_tm = mt.Tm_NN(seq, nn_table=mt.DNA_NN3, dnac1=50 if selfcomp else 25,
                      dnac2=0 if selfcomp else 25, selfcomp=selfcomp, Na=50, saltcorr=5)
    errors_bio.append(abs(bio_tm - point["tm"]))
    errors_p3.append(abs(primer3.calc_tm(seq, **conditions) - point["tm"]))
for point in data["picks"]:
    seq = point["oligo"]
    hairpin.append(primer3.calc_hairpin(seq, **conditions).tm)
    dimer.append(primer3.calc_homodimer(seq, **conditions).tm)
    assert abs(hairpin[-1] - point["hairpin"]) < 0.01
    assert abs(dimer[-1] - point["homodimer"]) < 0.01
    assert max(hairpin[-1], dimer[-1]) < 47
print(f"Biopython {Bio.__version__}, primer3-py {primer3.__version__}")
print(f"Independent oligos: {len(oligos)}; selected entries: {len(hairpin)}; plans: {data['counts']}")
print(f"Max absolute Tm error in C: Biopython {max(errors_bio):.8g}; Primer3 {max(errors_p3):.8g}")
print(f"Highest Primer3 structure Tm in selected entries: hairpin {max(hairpin):.2f} C, homodimer {max(dimer):.2f} C")
assert max(errors_bio) < 0.01 and max(errors_p3) < 0.01

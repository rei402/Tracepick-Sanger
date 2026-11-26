"use strict";
const test=require("node:test"), assert=require("node:assert/strict"), fs=require("node:fs"),path=require("node:path");
const S=require("../docs/engine.js");
const raw=fs.readFileSync(path.join(__dirname,"pUC19c-L09137.2.fa"),"utf8"), seq=S.parse(raw);
const bases=ranges=>ranges.flatMap(([a,b])=>Array.from({length:b-a},(_,i)=>a+i));
const normalize=(a,n)=>(a%n+n)%n;

test("single-record parser preserves coordinate numbering",()=>{
  assert.equal(S.parse(" >x\r\na cg\nt "),"ACGT");
  for(const bad of ["",">a\nA\n>b\nC","AC-N","AC1GT","ACU","A".repeat(20001)]) assert.throws(()=>S.parse(bad));
});
test("Tm and contiguous complementary structures",()=>{
  assert.equal(S.rc("ATCC"),"GGAT");
  assert.ok(S.selfComplementary("GCGCAAAAAGCGC"));
  assert.ok(S.selfComplementary("ATATCCCCCCCCATAT"));
  assert.ok(!S.selfComplementary("AGCAATAAACCAGCCAGCCGGAAG"));
});
test("invalid settings and coordinates",()=>{
  for(const opt of [{mode:"x"},{first:0},{last:1501},{first:200,last:100},{mode:"target",start:0},
    {mode:"target",start:500,stop:1},{mode:"target",stop:3.5},{mode:"target",stop:2687}])
    assert.throws(()=>S.design(seq,opt));
});
test("returns fewer when the template or intended read is unsuitable",()=>{
  for(const dna of ["ACGT","A".repeat(1000),seq.slice(0,55)]) assert.equal(S.design(dna).picked.length,0);
  assert.equal(S.design(seq,{mode:"target",start:100,stop:1000}).picked.length,0);
  assert.equal(S.design(seq,{direction:"F",mode:"target",start:1,stop:1}).picked.length,0);
});
test("duplicated copies have no qualifying internal binding site",()=>{
  const segment=seq.slice(0,240);
  const r=S.design(segment+segment,{circular:true});
  assert.equal(r.picked.length,0);assert.ok(r.counts.duplicate>0);
});
test("small supported interval returns a nonempty partial list",()=>{
  const p=S.design(seq.slice(0,115),{direction:"F",mode:"target",start:106,stop:106}).picked;
  assert.ok(p.length>0 && p.length<5);
});
test("single substitutions inside a unique 3′ tail are detected on both strands",()=>{
  const fragment=seq.slice(0,88);
  const mutate=s=>s.slice(0,42)+({A:"T",T:"A",C:"G",G:"C"}[s[42]])+s.slice(43);
  for(const copy of [mutate(fragment),S.rc(mutate(fragment))]){
    const input=fragment+seq.slice(610,637)+copy;
    const strands=[input,S.rc(input)];
    let rejected=0, missedByTailAlone=0;
    for(const s of strands) for(let i=0;i<s.length;i++) for(let len=18;len<=25 && i+len<=s.length;len++){
      const oligo=s.slice(i,i+len),gc=100*[...oligo].filter(b=>b==="G"||b==="C").length/len;
      if(gc<40||gc>60||/([ACGT])\1{3}|([ACGT]{2})\2{3}/.test(oligo)||
        S.tm(oligo)<55||S.tm(oligo)>65||S.selfComplementary(oligo))continue;
      let tails=0,near=0;
      for(const t of strands)for(let j=0;j<t.length;j++){
        if(t.slice(j,j+12)===oligo.slice(-12))tails++;
        if(j+len<=t.length && [...oligo].filter((b,k)=>b!==t[j+k]).length<=1)near++;
      }
      if(tails>1 || near>1)rejected++;
      if(tails===1 && near>1)missedByTailAlone++;
    }
    assert.ok(missedByTailAlone>0);assert.equal(S.design(input).counts.duplicate,rejected);
  }
});
test("Primer3 rejects stable structures and fills vacancies",()=>{
  const r=S.design(seq,{mode:"target",start:600,stop:660});
  assert.ok(r.rejected>0);assert.equal(r.picked.length,5);
  assert.ok(r.picked.every(p=>p.hairpin<47 && p.homodimer<47));
});
for (const [title,opt] of [
  ["linear template",{}], ["circular template",{circular:true}],["one target",{mode:"target",start:680,stop:729}],
  ["cross-origin target",{circular:true,mode:"target",start:2675,stop:19}],
  ["reverse only",{direction:"R"}],["forward only",{direction:"F"}],
  ["single base read window",{first:60,last:60,mode:"target",start:720,stop:720}],
  ["short circle, one lap",{circular:true,last:1500,first:3}]
]) test(title+" reconstructs each ordering sequence, read and union",()=>{
  const input=title.startsWith("short") ? seq.slice(0,160) : seq;
  const n=input.length,r=S.design(input,opt),union=new Set(), strands={F:new Set(),R:new Set()};
  if(!title.startsWith("single base")) assert.equal(r.picked.length,5);
  else assert.ok(r.picked.length<=2);
  const templateList=[input,S.rc(input)];
  for(const p of r.picked) {
    const plus=bases(p.bind).map(i=>input[i]).join("");
    assert.equal(p.oligo,p.dir==="F"?plus:S.rc(plus));
    if(opt.direction) assert.equal(p.dir,opt.direction);
    const step=p.dir==="F"?1:-1, terminal=p.dir==="F"?bases(p.bind).at(-1):bases(p.bind)[0];
    const expected=[];
    for(let d=r.opt.first;d<=r.opt.last && d<=n-p.oligo.length;d++){
      const i=terminal+step*d;
      if(opt.circular || (i>=0 && i<n)) expected.push(normalize(i,n));
    }
    const actual=bases(p.read);
    assert.deepEqual(new Set(actual),new Set(expected));
    assert.equal(actual.length,expected.length);
    if(r.opt.mode==="target") assert.ok(bases(r.target).every(i=>actual.includes(i)));
    for(const x of actual){union.add(x);strands[p.dir].add(x);}
    let matches=0,tails=0;
    for(const original of templateList){
      const t=original+(opt.circular?original:"");
      for(let i=0;i<n;i++){
        if(t.slice(i,i+12)===p.oligo.slice(-12)) tails++;
        const window=t.slice(i,i+p.oligo.length);
        if(window.length===p.oligo.length && [...p.oligo].filter((x,j)=>x!==window[j]).length<=1) matches++;
      }
    }
    assert.equal(tails,1);assert.equal(matches,1);
  }
  assert.equal(r.covered,union.size);
  assert.equal(r.dual,[...strands.F].filter(x=>strands.R.has(x)).length);
  assert.deepEqual(new Set(bases(r.gaps)),new Set(Array.from({length:n},(_,i)=>i).filter(i=>!union.has(i))));
});

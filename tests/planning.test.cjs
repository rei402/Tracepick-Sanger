"use strict";
const test=require("node:test"),assert=require("node:assert/strict");
const fs=require("node:fs"),path=require("node:path");
const S=require("../docs/engine.js"),R=require("../docs/risks.js"),P=require("../docs/planning.js");
const ref=S.parse(fs.readFileSync(path.join(__dirname,"pUC19c-L09137.2.fa"),"utf8"));

test("explicit counts select alternatives beyond five and stop when exhausted",()=>{
  const options={mode:"target",start:720,stop:720,first:30,last:750};
  const many=S.design(ref,{...options,count:15});
  assert.equal(many.picked.length,15);
  assert.deepEqual(S.design(ref,{...options,count:1}).picked,[many.picked[0]]);
  for(const count of [0,-2,3.5,NaN,Infinity]) assert.throws(()=>S.design(ref,{count}));
  const scarce=S.design(ref.slice(0,115),{count:10000000,direction:"F"});
  assert.ok(scarce.picked.length>0 && scarce.picked.length<20);
});
test("global button settings can exceed five, reach closure, and expose real gaps",()=>{
  const opt={mode:"global",start:99999,stop:-1,count:0,first:30,last:750};
  const r=S.design(ref,{...opt,circular:true});
  assert.ok(r.picked.length>5); assert.equal(r.covered,ref.length);
  const seen=new Set();
  for(const p of r.picked){
    const before=seen.size;
    for(const [a,b] of p.read)for(let i=a;i<b;i++)seen.add(i);
    assert.ok(seen.size>before);
  }
  const oneWay=S.design(ref,{...opt,direction:"F"});
  assert.ok(oneWay.gaps.length>0);assert.ok(oneWay.picked.every(p=>p.dir==="F"));
  assert.equal(oneWay.gaps[0][0],0);
  assert.equal(S.design("ACGT".repeat(100),opt).picked.length,0);
});
test("global greedy union equals all thermo-passing candidate reads, including adjacent 3′ ends",()=>{
  const pool=Array.from({length:430},(_,i)=>({
    dir:i%2?"F":"R",end3:i%25,oligo:String(i),read:[[i%37,i%37+1]],penalty:(i%7)/2
  }));
  const thermo={calculate:id=>({hairpin:Number(id)%37===7?50:20,homodimer:0})};
  const out=P.select(pool,50,{mode:"global",circular:false},thermo);
  const actual=[...out.depth].flatMap((v,i)=>v?[i]:[]);
  assert.deepEqual(actual,Array.from({length:37},(_,i)=>i).filter(i=>i!==7));
});
test("risk masks agree with independent motif enumeration across both template topologies",()=>{
  const seq="GGGTGCTAGTAACGCAAAAAATCTCTCTCTCCGGGGCCACGTTAAAGCGTACGTGGGGCCCCGGG";
  for(const circular of [false,true]){
    const wanted=Object.fromEntries(Object.keys(R.labels).map(k=>[k,new Set()]));
    for(let start=0;start<seq.length;start++){
      for(const len of [6,10,15,16,17,18,19,20,21,22,23,24]){
        if(!circular && start+len>seq.length)continue;
        const bases=Array.from({length:len},(_,j)=>(start+j)%seq.length);
        const s=bases.map(i=>seq[i]).join("");
        const mark=key=>bases.forEach(i=>wanted[key].add(i));
        if(len===20 && (s.match(/[GC]/g)||[]).length>=14)mark("gc");
        if(len===6 && new Set(s).size===1)mark("run");
        if(len===10 && s[0]!==s[1] && /^(..)(?:\1){4}$/.test(s))mark("tandem");
        if(len>=15 && s.slice(0,6)===S.rc(s.slice(-6)))mark("hairpin");
      }
    }
    const got=R.scan(seq,circular,S.rc);
    for(const key of Object.keys(wanted)){
      const actual=new Set(got.filter(f=>f.type===key).flatMap(f=>Array.from({length:f.b-f.a},(_,i)=>f.a+i)));
      assert.deepEqual(actual,wanted[key]);
    }
  }
  assert.ok(R.scan("ACGTCGAAACGACGT",false,S.rc).some(f=>f.type==="hairpin" && f.a===0 && f.b===15));
});
test("each primer's risk annotations intersect the correct template coordinates",()=>{
  for(const circular of [false,true]){
    const r=S.design(ref,{mode:"global",circular,first:30,last:750}),n=ref.length;
    for(const p of r.picked){
      const bases=ranges=>ranges.flatMap(([a,b])=>Array.from({length:b-a},(_,i)=>a+i));
      const early=Array.from({length:29},(_,j)=>p.end3+(p.dir==="F"?1:-1)*(1+j));
      const zones={bind:bases(p.bind),read:bases(p.read),lead:early.filter(i=>circular||(i>=0&&i<n)).map(i=>(i+n)%n)};
      for(const [zone,coords] of Object.entries(zones))for(const type of Object.keys(R.labels)){
        const expected=new Set(coords.filter(i=>r.features.some(f=>f.type===type && i>=f.a && i<f.b)));
        const actual=new Set(bases(p.warnings.filter(w=>w.zone===zone && w.type===type).map(w=>[w.a,w.b])));
        assert.deepEqual(actual,expected);
      }
    }
  }
});

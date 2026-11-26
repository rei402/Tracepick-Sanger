"use strict";
const TracepickRisks = (() => {
  const labels={gc:"高 GC", hairpin:"疑似发卡", run:"同聚碱基", tandem:"二核苷酸重复"};
  const zones={bind:"结合区", lead:"起读前延伸区", read:"预计读段"};
  // Each track is the union of matching motif footprints on the entered strand.
  function scan(seq,circular,rc) {
    const n=seq.length, s=circular?seq+seq:seq;
    const diff=Object.fromEntries(Object.keys(labels).map(k=>[k,new Int32Array(n+1)]));
    function mark(type,start,len) {
      if(len>n || start+len>s.length) return;
      const d=diff[type],end=start+len;
      d[start]++; d[Math.min(n,end)]--;
      if(end>n) {d[0]++; d[end-n]--;}
    }
    let gc=0;
    for(let i=0;i<Math.min(19,s.length);i++) gc+=Number("GC".includes(s[i]));
    for(let i=0;i<n;i++) {
      if(i+20<=s.length && n>=20) {
        gc+=Number("GC".includes(s[i+19]));
        if(gc>=14) mark("gc",i,20);
        gc-=Number("GC".includes(s[i]));
      }
      if(s.slice(i,i+6)===s[i].repeat(6)) mark("run",i,6);
      const pair=s.slice(i,i+2);
      if(pair.length===2 && pair[0]!==pair[1] && s.slice(i,i+10)===pair.repeat(5)) mark("tandem",i,10);
      const arm=rc(s.slice(i,i+6));
      for(let loop=3;loop<=12;loop++) {
        const len=12+loop;
        if(len>n || i+len>s.length) break;
        if(s.startsWith(arm,i+6+loop)) mark("hairpin",i,len);
      }
    }
    const features=[];
    for(const [type,d] of Object.entries(diff)) {
      let active=0,start=0;
      for(let i=0;i<=n;i++) {
        const prev=active;active+=d[i];
        if(!prev && active) start=i;
        if(prev && !active) features.push({type,a:start,b:i});
      }
    }
    return features;
  }
  function annotate(features,regions) {
    const warnings=[];
    for(const [zone,parts] of Object.entries(regions)) {
      for(const [left,right] of parts) for(const f of features) {
        const a=Math.max(f.a,left), b=Math.min(f.b,right);
        if(a<b) warnings.push({zone,type:f.type,a,b});
      }
    }
    return warnings;
  }
  return {scan,annotate,labels,zones};
})();
if(typeof module!=="undefined") module.exports=TracepickRisks;

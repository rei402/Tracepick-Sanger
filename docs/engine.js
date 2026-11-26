/* All internal intervals are zero-based and half-open on the entered strand.
   A candidate's "end3" is the index of its terminal nucleotide. */
"use strict";
const Sanger = (() => {
  const comp = { A: "T", C: "G", G: "C", T: "A" };
  const mod = (x, n) => ((x % n) + n) % n;
  const rc = seq => Array.from(seq, base => comp[base]).reverse().join("");
  const gcCount = seq => (seq.match(/[GC]/g) || []).length;
  const filters = {
    composition: "GC、同聚碱基或串联重复",
    tm: "Tm 超出范围",
    structure: "自身互补",
    duplicate: "模板上有其他可能结合位点",
    geometry: "读段位置不满足目标"
  };
  function parse(raw) {
    const lines = String(raw).trim().split(/\r\n?|\n/);
    if (lines[0].startsWith(">")) lines.shift();
    if (lines.some(line => line.trimStart().startsWith(">"))) throw Error("每次请输入一条 FASTA 记录。");
    const seq = lines.join("").replace(/\s/g, "").toUpperCase();
    if (!seq.length) throw Error("请粘贴 DNA 序列或导入 FASTA 文件。");
    const invalid = /[^ACGT]/.exec(seq);
    if (invalid) throw Error("第 " + (invalid.index + 1) + " 位为 " + invalid[0] + "。请提供仅含 A、C、G、T 的已确定模板。");
    if (seq.length > 20000) throw Error("当前支持最长 20,000 bp 的单条模板。");
    return seq;
  }

  // SantaLucia 1998 Table 2; delta H (kcal/mol), delta S (cal/mol/K).
  const nn = {};
  for (const [keys, h, s] of [
    ["AA TT",-7.9,-22.2],["AT",-7.2,-20.4],["TA",-7.2,-21.3],
    ["CA TG",-8.5,-22.7],["GT AC",-8.4,-22.4],["CT AG",-7.8,-21.0],
    ["GA TC",-8.2,-22.2],["CG",-10.6,-27.2],["GC",-9.8,-24.4],["GG CC",-8.0,-19.9]
  ]) for (const k of keys.split(" ")) nn[k] = [h, s];
  function tm(seq) {
    let h = 0, s = 0;
    for (const base of [seq[0], seq.at(-1)]) {
      const at = base === "A" || base === "T";
      h += at ? 2.3 : 0.1;
      s += at ? 4.1 : -2.8;
    }
    for (let i = 1; i < seq.length; i++) {
      const [dh, ds] = nn[seq.slice(i-1, i+1)];
      h += dh; s += ds;
    }
    s += 0.368 * (seq.length-1) * Math.log(0.05);
    const symmetric = seq === rc(seq);
    if (symmetric) s -= 1.4;
    const effective = symmetric ? 50e-9 : 50e-9 / 4;
    return 1000*h / (s + 1.987*Math.log(effective)) - 273.15;
  }
  function selfComplementary(seq) {
    const rev = rc(seq);
    // Terminal homodimer: at least four contiguous base pairs at the 3' end.
    if (seq.includes(rev.slice(0,4))) return true;
    for (let i = 0; i <= seq.length-6; i++)
      if (rev.includes(seq.slice(i,i+6))) return true;
    // Four paired bases on each arm and an intervening loop of >= 3 bases.
    for (let i = 0; i <= seq.length-11; i++)
      if (seq.indexOf(rc(seq.slice(i,i+4)), i+7) >= 0) return true;
    return false;
  }
  function parts(start, length, n, circular) {
    if (length <= 0) return [];
    if (!circular) {
      const left = Math.max(0,start), right = Math.min(n,start+length);
      return right > left ? [[left,right]] : [];
    }
    const left = mod(start,n), right = left+length;
    return right <= n ? [[left,right]] : [[left,n],[0,right-n]];
  }
  const width = ranges => ranges.reduce((sum,[a,b]) => sum+b-a,0);
  function intersection(first, second) {
    let bp = 0;
    for (const [a,b] of first) for (const [c,d] of second)
      bp += Math.max(0, Math.min(b,d)-Math.max(a,c));
    return bp;
  }
  function reference(seq, circular) {
    const n = seq.length;
    const orientations = [seq,rc(seq)].map(s => circular ? s+s : s);
    const seeds = new Map(), tails = new Map();
    for (let strand = 0; strand < 2; strand++) {
      const template = orientations[strand];
      for (let i=0; i<n; i++) {
        if (i+8<=template.length) {
          const key = template.slice(i,i+8);
          if (!seeds.has(key)) seeds.set(key,[]);
          seeds.get(key).push([strand,i]);
        }
        if (i+12<=template.length) {
          const tail = template.slice(i,i+12);
          tails.set(tail,(tails.get(tail)||0)+1);
        }
      }
    }
    return { n, circular, orientations, seeds, tails };
  }
  function unique(oligo, strand, pos, ref) {
    if (ref.tails.get(oligo.slice(-12)) !== 1) return false;
    // For <=1 substitution at least one of two disjoint 8-base seeds is exact.
    for (const offset of [0, oligo.length-8]) {
      const key = oligo.slice(offset,offset+8);
      for (const [other,hit] of ref.seeds.get(key) || []) {
        let start = hit-offset;
        if (ref.circular) start=mod(start,ref.n);
        if (other===strand && start===pos) continue;
        const s = ref.orientations[other];
        if (start<0 || start+oligo.length>s.length) continue;
        let mismatches = 0;
        for (let i=0; i<oligo.length && mismatches<2; i++)
          mismatches += Number(oligo[i]!==s[start+i]);
        if (mismatches<=1) return false;
      }
    }
    return true;
  }

  function design(raw, settings={}) {
    const thermo = typeof SangerThermo!=="undefined" ? SangerThermo : require("./thermo.js");
    const seq=parse(raw), n=seq.length;
    const opt={circular:false, direction:"both", mode:"coverage", start:1, stop:n,
      first:50, last:500, ...settings};
    if (!["both","F","R"].includes(opt.direction) || !["coverage","target"].includes(opt.mode) ||
        typeof opt.circular !== "boolean") throw Error("未知的设计设置。");
    if (![opt.first,opt.last].every(Number.isSafeInteger) || opt.first<1 || opt.last<opt.first || opt.last>1500)
      throw Error("预计读段的起止距离须为 1–1,500 的整数，且起点不大于终点。");
    let target = [[0,n]];
    if (opt.mode==="target") {
      if (![opt.start,opt.stop].every(x=>Number.isSafeInteger(x)&&x>=1&&x<=n) ||
          (!opt.circular && opt.start>opt.stop)) throw Error("请检查目标坐标，范围为 1–"+n+"；只有环状模板可以跨越原点。");
      target=opt.start<=opt.stop ? [[opt.start-1,opt.stop]] : [[opt.start-1,n],[0,opt.stop]];
    }
    const ref=reference(seq,opt.circular), counts=Object.fromEntries(Object.keys(filters).map(k=>[k,0]));
    const pool=[];
    for (let strand=0; strand<2; strand++) {
      const dir=strand===0 ? "F" : "R";
      if (opt.direction!=="both" && dir!==opt.direction) continue;
      const oriented=ref.orientations[strand];
      for (let pos=0;pos<n;pos++) for (let len=18;len<=Math.min(25,n);len++) {
        if (!opt.circular && pos+len>n) break;
        const oligo=oriented.slice(pos,pos+len);
        const gc=100*gcCount(oligo)/len;
        if (gc<40 || gc>60 || /([ACGT])\1{3}|([ACGT]{2})\2{3}/.test(oligo)) {counts.composition++;continue;}
        const temp=tm(oligo);
        if (temp<55 || temp>65) {counts.tm++;continue;}
        if (selfComplementary(oligo)) {counts.structure++;continue;}
        if (!unique(oligo,strand,pos,ref)) {counts.duplicate++;continue;}
        const left = strand===0 ? pos : n-pos-len;
        const end3 = mod(strand===0 ? left+len-1 : left,n);
        const far = Math.min(opt.last,n-len);
        const read = parts(strand===0 ? end3+opt.first : end3-far,
          far-opt.first+1,n,opt.circular);
        if (!read.length || (opt.mode==="target" && intersection(read,target)!==width(target))) {
          counts.geometry++; continue;
        }
        const penalty = Math.abs(temp-60) + 0.08*Math.abs(gc-50) + 0.18*Math.abs(len-21)
          + (/[GC]$/.test(oligo) ? 0 : 0.4);
        pool.push({dir,end3,oligo,tm:temp,gc,penalty,bind:parts(left,len,n,opt.circular),read});
      }
    }
    const picked=[], depth=new Uint8Array(n), byDir={F:new Uint8Array(n),R:new Uint8Array(n)};
    const rejected=new Set(), tested=new Map();
    while (picked.length<5) {
      function uncoveredPrefix(arr) {
        const prefix = new Uint32Array(n+1);
        for (let i=0;i<n;i++) prefix[i+1]=prefix[i]+Number(arr[i]===0);
        return prefix;
      }
      const all=uncoveredPrefix(depth), forward=uncoveredPrefix(byDir.F), reverse=uncoveredPrefix(byDir.R);
      const ranked=[];
      for (const p of pool) {
        if (rejected.has(p.oligo)) continue;
        if (picked.some(q => {
          const dist=Math.abs(q.end3-p.end3);
          return p.dir===q.dir && (opt.circular ? Math.min(dist,n-dist) : dist)<12;
        })) continue;
        let fresh=0, strandFresh=0;
        const same=p.dir==="F" ? forward : reverse;
        for (const [a,b] of p.read) {
          fresh+=all[b]-all[a]; strandFresh+=same[b]-same[a];
        }
        // Lexicographic priorities; bounded penalty is < 10 and read length <= 1500.
        const utility=opt.mode==="coverage" ? fresh*20000+strandFresh*10-p.penalty :
          (picked.some(q=>q.dir===p.dir) ? 0 : 20)-p.penalty;
        ranked.push({p,utility});
      }
      ranked.sort((a,b)=>b.utility-a.utility);
      let best=null;
      for(const {p} of ranked) {
        if(!tested.has(p.oligo)) tested.set(p.oligo,thermo.calculate(p.oligo));
        const structure=tested.get(p.oligo);
        if(structure.hairpin>=47 || structure.homodimer>=47) {rejected.add(p.oligo);continue;}
        best={...p,...structure};break;
      }
      if (!best) break;
      picked.push(best);
      for (const [a,b] of best.read) for (let i=a;i<b;i++) {depth[i]++;byDir[best.dir][i]++;}
    }
    const gaps=[];
    for (let i=0;i<n;) {
      if (depth[i]) {i++;continue;}
      const start=i;
      while(i<n && !depth[i]) i++;
      gaps.push([start,i]);
    }
    return {n,opt,target,picked,counts,pass:pool.length,gaps,checked:tested.size,rejected:rejected.size,covered:n-width(gaps),
      dual:byDir.F.reduce((total,val,i)=>total+Number(val>0 && byDir.R[i]>0),0)};
  }
  return {parse,tm,rc,design,filters,selfComplementary};
})();
if (typeof module!=="undefined" && module.exports) module.exports=Sanger;

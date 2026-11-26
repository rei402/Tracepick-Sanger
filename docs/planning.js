"use strict";
const TracepickPlanning = (() => {
  // Coverage grows monotonically. Each cached candidate score is an upper bound.
  // Break ties in enumeration order, including after a stale score is recomputed.
  function better(a,b) {
    return a.score>b.score || (a.score===b.score && a.order<b.order);
  }
  function siftDown(heap,index) {
    const item=heap[index];
    let child=2*index+1;
    while(child<heap.length) {
      if(child+1<heap.length && better(heap[child+1],heap[child])) child++;
      if(!better(heap[child],item)) break;
      heap[index]=heap[child]; index=child; child=2*index+1;
    }
    heap[index]=item;
  }
  function coverage(n) {
    const seen=new Uint8Array(n), tree=new Int32Array(n+1);
    function prefix(i) {
      let sum=0;
      while(i>0) {sum+=tree[i]; i-=i&-i;}
      return sum;
    }
    return {
      seen,
      fresh(ranges) {
        let count=0;
        for(const [a,b] of ranges) count+=b-a-prefix(b)+prefix(a);
        return count;
      },
      mark(ranges) {
        for(const [a,b] of ranges) for(let i=a;i<b;i++) {
          if(seen[i]) continue;
          seen[i]=1;
          for(let j=i+1;j<=n;j+=j&-j) tree[j]++;
        }
      }
    };
  }
  function select(pool,n,opt,thermo) {
    const all=coverage(n), strands={F:coverage(n),R:coverage(n)};
    const blocked={F:new Uint8Array(n),R:new Uint8Array(n)}, usedDir=new Set();
    const picked=[];
    let checked=0, rejected=0;
    function score(p) {
      if(opt.mode==="target") return (usedDir.has(p.dir)?0:20)-p.penalty;
      const fresh=all.fresh(p.read);
      if(opt.mode==="global" && fresh===0) return -Infinity;
      return fresh*20000+strands[p.dir].fresh(p.read)*10-p.penalty;
    }
    const heap=pool.map((p,order)=>({p,order,score:score(p)}));
    for(let i=Math.floor(heap.length/2)-1;i>=0;i--) siftDown(heap,i);
    while(heap.length && (opt.mode==="global" || picked.length<opt.count)) {
      const top=heap[0],p=top.p;
      if(!blocked[p.dir][p.end3]) {
        const current=score(p);
        if(current!==top.score) {
          top.score=current; siftDown(heap,0); continue;
        }
        if(current===-Infinity) break;
        const structure=thermo.calculate(p.oligo);
        checked++;
        if(structure.hairpin<47 && structure.homodimer<47) {
          picked.push({...p,...structure});usedDir.add(p.dir);
          all.mark(p.read);strands[p.dir].mark(p.read);
          if(opt.mode!=="global") for(let d=-11;d<=11;d++) {
            const i=p.end3+d;
            if(opt.circular) blocked[p.dir][(i%n+n)%n]=1;
            else if(i>=0 && i<n) blocked[p.dir][i]=1;
          }
        } else rejected++;
      }
      const last=heap.pop();
      if(heap.length) {heap[0]=last;siftDown(heap,0);}
    }
    return {picked,checked,rejected,depth:all.seen,byDir:{F:strands.F.seen,R:strands.R.seen}};
  }
  return {select};
})();
if(typeof module!=="undefined") module.exports=TracepickPlanning;

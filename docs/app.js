"use strict";
(() => {
  const el = id => document.getElementById(id);
  const targetHelp = el("target-note").textContent;
  const formatRanges = ranges => ranges.map(([a,b]) => a+1===b ? String(b) : (a+1)+"–"+b).join(" / ");
  const named = (p,i) => p.dir+String(i+1).padStart(2,"0");
  const percent = (a,b) => (100*a/b).toFixed(1)+"%";
  let report=null, revision=0;
  function clear() {
    revision++;
    report=null; el("result").hidden=true; el("status").textContent="";
    const wholeTemplate=el("mode").value==="coverage";
    el("start").disabled=wholeTemplate;
    el("stop").disabled=wholeTemplate;
    el("target-note").textContent=wholeTemplate ? "当前按整段模板设计。切回指定位点模式可继续使用这组坐标。" : targetHelp;
    try {el("length").textContent=Sanger.parse(el("dna").value).length.toLocaleString()+" bp";}
    catch {el("length").textContent=el("dna").value.trim() ? "待检查" : "0 bp";}
  }
  el("form").addEventListener("input", clear);
  el("example").onclick=()=>{
    // Reproducible artificial sequence for the interface demonstration.
    let state=70211, seq="";
    for(let i=0;i<1450;i++){
      state=(Math.imul(state,1664525)+1013904223)>>>0;
      seq+="ACGT"[state>>>30];
    }
    el("dna").value=">demo_1450bp\n"+seq.match(/.{1,70}/g).join("\n");
    el("start").value=720; el("stop").value=720;
    clear();
  };
  el("import").onclick=()=>el("file").click();
  el("file").onchange=async()=>{
    const file=el("file").files[0];
    if (!file) return;
    clear(); const before=revision;
    try {
      if(file.size>300000) throw Error("文件过大。请导入一条 20,000 bp 以内的 FASTA 序列。");
      const value=await file.text();
      if (revision!==before) return;
      Sanger.parse(value); el("dna").value=value; clear();
    } catch(e) { if (before===revision) el("status").textContent=e.message; }
    finally {el("file").value="";}
  };
  function element(tag, text, cls) {
    const node=document.createElement(tag);
    if (text!==undefined) node.textContent=text;
    if (cls) node.className=cls;
    return node;
  }
  function lane(label, read, bind, type, n) {
    const row=element("div",undefined,"lane"), track=element("div",undefined,"track");
    row.append(element("span",label),track);
    for(const [ranges,cls] of [[read,type],[bind,"bind"]]) for(const [a,b] of ranges){
      const bar=element("span",undefined,"segment "+cls);
      bar.style.left=(100*a/n)+"%"; bar.style.width=(100*(b-a)/n)+"%";
      bar.title=label+" · "+formatRanges([[a,b]]); track.append(bar);
    }
    return row;
  }
  function render(r) {
    const count=r.picked.length;
    el("num").textContent=count+" / 5";
    el("summary").textContent=r.opt.mode==="coverage" ?
      "预计合计覆盖 "+r.covered.toLocaleString()+" / "+r.n.toLocaleString()+" bp（"+percent(r.covered,r.n)+"），其中 "+r.dual+" bp 可从两个方向读取。" :
      "目标 "+formatRanges(r.target)+" bp · "+count+" 条候选可分别读完整个目标。";
    const map=el("map"), head=element("div",undefined,"maphead");
    head.append(element("strong","预计读段覆盖"),element("span","深色为引物 · 绿色正向 → · 橙色反向 ←","legend"));
    const scale=element("div",undefined,"scale"),axis=element("div",undefined,"axis");
    axis.append(element("span","1 bp"),element("span",r.n.toLocaleString()+" bp"));scale.append(element("span",""),axis);
    map.replaceChildren(head,scale);
    if(r.opt.mode==="target") map.append(lane("目标",r.target,[],"target",r.n));
    const rows=el("rows"); rows.replaceChildren();
    r.picked.forEach((p,i)=>{
      const name=named(p,i), row=element("tr"), cell=element("td"),button=element("button",p.oligo,"seq");
      button.type="button";button.title="复制 "+name+" 的 5′ → 3′ 序列";
      button.setAttribute("aria-label",button.title+" "+p.oligo);
      button.onclick=async()=>{
        try {await navigator.clipboard.writeText(p.oligo);el("status").textContent="已复制 "+name+" 的订购序列。";}
        catch {const range=document.createRange();range.selectNodeContents(button);
          const selection=window.getSelection();selection.removeAllRanges();selection.addRange(range);
          el("status").textContent="已选中序列，请使用系统复制命令。";}
      };
      cell.append(button);
      row.append(element("td",name+(p.dir==="F"?" →":" ←")),cell);
      for(const value of [formatRanges(p.bind),p.oligo.length,p.tm.toFixed(2),p.gc.toFixed(1),formatRanges(p.read)])
        row.append(element("td",value));
      rows.append(row);map.append(lane(name,p.read,p.bind,p.dir,r.n));
    });
    const stats=el("stats");stats.replaceChildren();
    const evaluated=Object.values(r.counts).reduce((a,b)=>a+b,r.pass);
    stats.append(element("p","共检查 "+evaluated.toLocaleString()+" 个定向窗口，"+r.pass.toLocaleString()+" 个进入排序。下列初筛统计按首个排除原因计数。"));
    stats.append(element("p",Object.entries(r.counts).map(([k,v])=>Sanger.filters[k]+" "+v.toLocaleString()).join(" · ")));
    stats.append(element("p","按排名使用 Primer3 复核 "+r.checked+" 条序列，其中 "+r.rejected+" 条因发卡或同源二聚体 Tm ≥ 47 °C 被排除。"));
    for (const [i,p] of r.picked.entries()) stats.append(element("p",named(p,i)+"　发卡 Tm "+p.hairpin.toFixed(2)+" °C　·　同源二聚体 Tm "+p.homodimer.toFixed(2)+" °C"));
    stats.append(element("p","模板未覆盖区（bp）　"+(formatRanges(r.gaps)||"无")));
    stats.append(element("p","坐标为输入链上的 1-based 闭区间；斜线表示跨越环状原点。反向引物已转为可订购的 5′ → 3′ 序列。"));
    if(count<5){
      let reason="";
      const sorted=Object.entries(r.counts).sort((a,b)=>b[1]-a[1]);
      if(r.n<18) reason="模板长度不足最短引物长度 18 nt。";
      else if(r.pass) reason="经 Primer3 结构复核，并要求同方向引物的 3′ 端至少相距 12 bp，最终可选出 "+count+" 条。详见筛选统计。";
      else if(r.opt.mode==="target" && r.target.reduce((s,[a,b])=>s+b-a,0)>r.opt.last-r.opt.first+1)
        reason="目标长于一次反应的预计可读窗口。可改用整段覆盖模式。";
      else if(r.n-18<r.opt.first) reason="模板扣除最短引物后，剩余序列不足以达到设定的读段起始距离。";
      else reason="没有通过全部条件的候选。排除最多的条件为“"+Sanger.filters[sorted[0][0]]+"”。";
      el("status").textContent=reason;
    }
    el("export").disabled=!count;el("result").hidden=false;
  }
  el("form").onsubmit=async event=>{
    event.preventDefault();clear(); const before=revision;
    el("go").disabled=true;el("status").textContent="正在筛选模板上的候选位点…";
    await new Promise(resolve=>setTimeout(resolve,25));
    try {
      if(revision!==before) return;
      report=Sanger.design(el("dna").value,{
        mode:el("mode").value,circular:el("topology").value==="circular",direction:el("dir").value,
        start:el("start").valueAsNumber,stop:el("stop").valueAsNumber,
        first:el("first").valueAsNumber,last:el("last").valueAsNumber
      });
      el("status").textContent="";render(report);
      el("result").scrollIntoView({block:"nearest"});
    } catch(e) {el("status").textContent=e.message;}
    finally {el("go").disabled=false;}
  };
  el("export").onclick=()=>{
    if(!report) return;
    const r=report;
    const header=["name","sequence_5to3","direction","length_nt","Tm_C","GC_pct","binding_1based","read_1based","template_bp","topology","mode","target_1based","first_offset","last_offset","hairpin_Tm_C","homodimer_Tm_C","designer_version"];
    const lines=[header,...r.picked.map((p,i)=>[named(p,i),p.oligo,p.dir,p.oligo.length,p.tm.toFixed(2),p.gc.toFixed(1),
      formatRanges(p.bind),formatRanges(p.read),r.n,r.opt.circular?"circular":"linear",r.opt.mode,
      formatRanges(r.target),r.opt.first,r.opt.last,p.hairpin.toFixed(2),p.homodimer.toFixed(2),"1.1"])];
    const url=URL.createObjectURL(new Blob(["\uFEFF"+lines.map(x=>x.join("\t")).join("\r\n")],{type:"text/tab-separated-values;charset=utf-8"}));
    const link=element("a"); link.href=url;link.download="sanger-primers.tsv";document.body.append(link);link.click();link.remove();
    setTimeout(()=>URL.revokeObjectURL(url),1000);
  };
  clear();
})();

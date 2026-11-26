"use strict";
const SangerThermo = (() => {
  let api;
  function initialize() {
    const data=typeof THAL_WASM==="string" ? THAL_WASM : require("./thermo-data.js");
    const bytes=Uint8Array.from(atob(data),c=>c.charCodeAt(0));
    const mod=new WebAssembly.Module(bytes);
    let instance;
    const wasi={
      fd_close:()=>8, fd_seek:()=>8, fd_write:()=>8,
      fd_fdstat_get:()=>8, fd_prestat_get:()=>8, fd_prestat_dir_name:()=>8,
      random_get:(ptr,len)=>{
        globalThis.crypto.getRandomValues(new Uint8Array(instance.exports.memory.buffer,ptr,len));
        return 0;
      },
      proc_exit:()=>{throw Error("Primer3 计算异常终止。");}
    };
    instance=new WebAssembly.Instance(mod,{wasi_snapshot_preview1:wasi});
    instance.exports._initialize();
    if(instance.exports.prepare()!==0) throw Error("Primer3 参数初始化失败。");
    return instance.exports;
  }
  function calculate(seq) {
    if(!/^[ACGT]{1,25}$/.test(seq)) throw Error("结构计算需要 1–25 nt 的确定 DNA 序列。");
    try {
      api ||= initialize();
      const view=new Uint8Array(api.memory.buffer,api.input_address(),26);
      for(let i=0;i<seq.length;i++) view[i]=seq.charCodeAt(i);
      view[seq.length]=0;
      const hairpin=api.structure_tm(1), homodimer=api.structure_tm(0);
      if(!Number.isFinite(hairpin) || !Number.isFinite(homodimer)) throw Error("引擎返回异常数值。");
      return {hairpin,homodimer};
    } catch(error) {
      api=undefined;
      throw Error("Primer3 结构计算失败，请使用支持 WebAssembly 的现代浏览器。"+error.message);
    }
  }
  return {calculate};
})();
if(typeof module!=="undefined") module.exports=SangerThermo;

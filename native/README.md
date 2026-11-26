# Primer3 thermodynamics module

`upstream/` contains the THAL C implementation and embedded parameter tables distributed with primer3-py 2.3.0, plus the GPL text. The original notices are in each source file. The upstream files are included verbatim. Tracepick's C bridge sets the input buffer and applies oligo defaults (50 mM monovalent salt, 0 Mg and dNTP, 50 nM DNA). It invokes the hairpin and unrestricted homodimer calculations on the same oligo. [primer3-py source](https://github.com/libnano/primer3-py) · [Primer3 source](https://github.com/primer3-org/primer3)

The executable is an offline WebAssembly reactor built using Zig 0.13.0's C compiler and WASI libc. The small header in `shim/` converts fatal upstream longjmps into a WebAssembly trap. On failure the JavaScript interface discards the instance and reports an error. No file system or network access is provided to THAL.

## Rebuild

Use Zig 0.13.0 and Python 3. From the repository root run:

```sh
python native/build.py
```

Set the environment variable `ZIG` to the Zig executable if it is outside PATH. The script writes `docs/thermo-data.js`. The native-source zip on the site contains the build inputs and JavaScript interface with the same directory layout. The generated base64 representation permits loading the module from local files without a development server. The browser interface requires WebAssembly and a modern JavaScript runtime.

## Licenses

The application and THAL source are distributed under GPL-2.0-or-later. The build links portions of the Zig compiler runtime and WASI libc. Their notices are under `notices/`. [Zig 0.13.0 sources](https://github.com/ziglang/zig/tree/0.13.0)

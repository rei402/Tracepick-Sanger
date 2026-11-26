"""Build the bundled Primer3 module with Zig 0.13.0, Python 3 and no network."""
import base64
import os
import subprocess
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parent
with tempfile.TemporaryDirectory() as tmp:
    wasm = Path(tmp) / "thermo.wasm"
    subprocess.run([os.environ.get("ZIG", "zig"), "cc", "-target", "wasm32-wasi",
                    "-mexec-model=reactor", "-O2", "-ffunction-sections", "-fdata-sections",
                    "-Ishim", "-Iupstream", "bridge.c", "upstream/thal.c", "upstream/thal_parameters.c",
                    "-Wl,--export=input_address", "-Wl,--export=prepare", "-Wl,--export=structure_tm",
                    "-Wl,--strip-all", "-o", str(wasm)], cwd=root, check=True)
    encoded = base64.b64encode(wasm.read_bytes()).decode("ascii")
target = root.parent / "docs/thermo-data.js"
target.parent.mkdir(exist_ok=True)
target.write_text('/* Primer3 THAL. Source and build instructions: native/ */\n'
                  'const THAL_WASM="' + encoded + '";\n'
                  'if(typeof module!=="undefined") module.exports=THAL_WASM;\n', encoding="ascii")

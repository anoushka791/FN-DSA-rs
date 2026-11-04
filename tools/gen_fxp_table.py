# Regenerates src/fxp_tables.rs from ../kgen_fxp.c (run from the rust/ directory).
import re
src = open('../kgen_fxp.c').read()
m = re.search(r'static const fxc GM_TAB\[1024\] = \{(.*?)\n\};', src, re.S)
rows = re.findall(r'FXC\(\s*(\d+)ull,\s*(\d+)ull\)', m.group(1))
assert len(rows) == 1024, len(rows)
out = ["// Generated from kgen_fxp.c by tools/gen_fxp_table.py -- do not edit.\n\n",
       "/// Primitive 2048-th roots of unity (real, imaginary) in 32.32 fixed point,\n",
       "/// in the order needed by the FFT.\n",
       "pub static GM_TAB: [(u64, u64); 1024] = [\n"]
for i in range(0, 1024, 2):
    out.append("    " + ", ".join("(%sull, %sull)" % rows[j] for j in (i, i + 1)) + ",\n")
out.append("];\n")
s = "".join(out).replace("ull", "u64")
open('src/fxp_tables.rs', 'w').write(s)

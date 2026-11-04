# Regenerates src/primes.rs from ../kgen_mp31.c (run from the rust/ directory).
import re
src = open('../kgen_mp31.c').read()
m = re.search(r'const small_prime PRIMES\[\] = \{(.*?)\n\};', src, re.S)
rows = re.findall(r'\{\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+)\s*\}', m.group(1))
assert len(rows) >= 308, len(rows)
out = ["// Generated from kgen_mp31.c by tools/gen_primes.py -- do not edit.\n\n",
       "use crate::mp31::SmallPrime;\n\n",
       "pub static PRIMES: [SmallPrime; %d] = [\n" % len(rows)]
for r in rows:
    out.append("    SmallPrime { p: %s, p0i: %s, r2: %s, g: %s, ig: %s, s: %s },\n" % r)
out.append("];\n")
open('src/primes.rs', 'w').writelines(out)

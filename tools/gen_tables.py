# Regenerates src/tables.rs from ../mq.c (run from the rust/ directory).
import re
src=open('../mq.c').read()
def table(name):
    m=re.search(r'const uint16_t '+name+r'\[\] = \{(.*?)\};',src,re.S)
    return [int(x) for x in re.findall(r'\b\d+\b',m.group(1))]
out=["// Generated from mq.c by tools/gen_tables.py -- do not edit by hand.\n",
     "#[repr(C, align(32))]\npub struct Aligned32<T>(pub T);\n\n"]
for n,rn in (('mq_GM','GM'),('mq_iGM','IGM')):
    v=table(n); assert len(v)==1024
    out.append("pub static %s: Aligned32<[u16; 1024]> = Aligned32([\n"%rn)
    for i in range(0,1024,16):
        out.append("    "+", ".join(map(str,v[i:i+16]))+",\n")
    out.append("]);\n\n")
open('src/tables.rs','w').writelines(out)

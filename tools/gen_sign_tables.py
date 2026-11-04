# Regenerates src/sign_tables.rs from ../sign_fpoly.c and ../sign_mq2.c
# (run from the rust/ directory).
import re
import struct
import math


def fpr_bits(i, e):
    if i == 0:
        return 0
    v = math.ldexp(float(i), e)
    assert int(v / math.ldexp(1.0, e)) == i, (i, e)  # exactness check
    return struct.unpack('<Q', struct.pack('<d', v))[0]


def fpr_table(src, name):
    m = re.search(r'(?:static )?const fpr ' + name + r'\[\] = \{(.*?)\n\};', src, re.S)
    body = re.sub(r'/\*.*?\*/', '', m.group(1), flags=re.S)
    out = []
    for tok in re.finditer(r'FPR_ZERO|FPR\(\s*(-?\d+)\s*,\s*(-?\d+)\s*\)', body):
        if tok.group(0) == 'FPR_ZERO':
            out.append(0)
        else:
            out.append(fpr_bits(int(tok.group(1)), int(tok.group(2))))
    return out


def u16_table(src, name):
    m = re.search(r'const uint16_t ' + name + r'\[\] = \{(.*?)\};', src, re.S)
    return [int(x) for x in re.findall(r'\b\d+\b', m.group(1))]


fp = open('../sign_fpoly.c').read()
m2 = open('../sign_mq2.c').read()

out = ["// Generated from sign_fpoly.c / sign_mq2.c by tools/gen_sign_tables.py -- do not edit.\n\n",
       "#![allow(clippy::unreadable_literal)]\n\n",
       "use crate::tables::Aligned32;\n\n"]

for cname, rname, size in (('fpoly_SIN', 'SIN', 512), ('RNJ', 'RNJ', 512), ('IRNJ', 'IRNJ', 512),
                           ('X_PLUS_1_SQUARED_OVER_X', 'X_PLUS_1_SQUARED_OVER_X', 512),
                           ('X_MINUS_1_OVER_X_PLUS_1_TIMES_I', 'X_MINUS_1_OVER_X_PLUS_1_TIMES_I', 512)):
    t = fpr_table(fp, cname)
    assert len(t) == size, (cname, len(t))
    out.append("pub static %s: [f64; %d] = [\n" % (rname, size))
    for i in range(0, size, 4):
        out.append("    " + ", ".join("f64::from_bits(0x%016X)" % v for v in t[i:i + 4]) + ",\n")
    out.append("];\n\n")

for cname, rname in (('mq2_GM', 'GM2'), ('mq2_iGM', 'IGM2')):
    v = u16_table(m2, cname)
    assert len(v) == 1024, (cname, len(v))
    out.append("pub static %s: Aligned32<[u16; 1024]> = Aligned32([\n" % rname)
    for i in range(0, 1024, 16):
        out.append("    " + ", ".join(map(str, v[i:i + 16])) + ",\n")
    out.append("]);\n\n")

open('src/sign_tables.rs', 'w').writelines(out)

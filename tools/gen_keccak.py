# Mechanically translates the Keccak-f[1600] round body from ../sha3.c
# into src/keccak_gen.rs (run from the rust/ directory).
import re
lines=open('../sha3.c').read().split('\n')
# loop body: between 'for (j = 0; j < 24; j += 2) {' and its closing brace
start=next(i for i,l in enumerate(lines) if 'for (j = 0; j < 24; j += 2)' in l)+1
end=next(i for i in range(start,len(lines)) if lines[i]=='\t}')
body=lines[start:end]
out=[]
for l in body:
    l=l.strip()
    if not l: out.append(''); continue
    l=re.sub(r'A\[\s*(\d+)\]',r'a[\1]',l)
    l=re.sub(r'RC\[j \+ (\d+)\]',r'RC[j + \1]',l)
    l=l.replace('~','!')
    assert l.endswith(';'),l
    out.append('        '+l)
hdr='''// Generated from sha3.c by tools/gen_keccak.py -- do not edit by hand.
// Keccak-f[1600] with the "lane complementing" internal representation
// (words 1, 2, 8, 12, 17, 20 are stored inverted), two rounds per iteration.

const RC: [u64; 24] = [
    0x0000000000000001, 0x0000000000008082, 0x800000000000808A, 0x8000000080008000,
    0x000000000000808B, 0x0000000080000001, 0x8000000080008081, 0x8000000000008009,
    0x000000000000008A, 0x0000000000000088, 0x0000000080008009, 0x000000008000000A,
    0x000000008000808B, 0x800000000000008B, 0x8000000000008089, 0x8000000000008003,
    0x8000000000008002, 0x8000000000000080, 0x000000000000800A, 0x800000008000000A,
    0x8000000080008081, 0x8000000000008080, 0x0000000080000001, 0x8000000080008008,
];

#[inline(never)]
pub fn process_block(a: &mut [u64; 25]) {
    let (mut t0, mut t1, mut t2, mut t3, mut t4): (u64, u64, u64, u64, u64);
    let (mut tt0, mut tt1, mut tt2, mut tt3): (u64, u64, u64, u64);
    let (mut t, mut kt): (u64, u64);
    let (mut c0, mut c1, mut c2, mut c3, mut c4, mut bnn): (u64, u64, u64, u64, u64, u64);

    a[1] = !a[1];
    a[2] = !a[2];
    a[8] = !a[8];
    a[12] = !a[12];
    a[17] = !a[17];
    a[20] = !a[20];

    let mut j = 0usize;
    while j < 24 {
'''
ftr='''        j += 2;
    }

    a[1] = !a[1];
    a[2] = !a[2];
    a[8] = !a[8];
    a[12] = !a[12];
    a[17] = !a[17];
    a[20] = !a[20];
}
'''
open('src/keccak_gen.rs','w').write(hdr+'\n'.join(out)+'\n'+ftr)

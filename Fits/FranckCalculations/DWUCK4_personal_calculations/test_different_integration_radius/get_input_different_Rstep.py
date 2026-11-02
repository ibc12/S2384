#!/usr/bin/env python3
# make_dr.py  entrada  salida  DRF_nuevo  RMAX
import sys, re
fin, fout = sys.argv[1], sys.argv[2]
drf_new, rmax = float(sys.argv[3]), float(sys.argv[4])
drf_old = 0.1

step = round(drf_new / drf_old)
if abs(drf_new / drf_old - step) > 1e-9 or step < 1:
    sys.exit("DRF nuevo debe ser multiplo entero de 0.1")
npt = round(rmax / drf_new)
if npt > 400 or npt % 5:
    sys.exit("RMAX/DRF=%d: debe ser <=400 y multiplo de 5" % npt)

lines = open(fin).read().rstrip('\n').split('\n')

# tarjeta 4 (4a linea): DRF, RZ, RMAX, VCE, FNRNG, AMASS en campos de 8
f = [lines[3].ljust(80)[i:i+8] for i in range(0, 80, 8)]
assert abs(float(f[0]) - drf_old) < 1e-9, "la tarjeta 4 no es la esperada"
f[0] = "%8.4f" % drf_new
f[2] = "%8.4f" % rmax
lines[3] = ''.join(f).rstrip()

hdr = re.compile(r'^\s*(\d+)\.\s+([0-3])\.\s*$')
out, i = [], 0
while i < len(lines):
    m = hdr.match(lines[i]) if i > 3 else None
    if not m:
        out.append(lines[i]); i += 1; continue
    n, comp = int(m.group(1)), int(m.group(2))
    nl = (n + 4) // 5
    toks = [t for l in lines[i+1:i+1+nl] for t in l.split()]
    nxt = lines[i+1+nl].split() if i+1+nl < len(lines) else []
    if len(toks) != n or (nxt and all('E' in t for t in nxt)):
        sys.exit("Tabla en linea %d: no tiene %d valores" % (i+1, n))
    last = max((j for j, t in enumerate(toks, 1) if abs(float(t)) > 1e-6), default=0)
    keep = toks[step-1::step][:npt]
    keep += ['0.0000000E+00'] * (npt - len(keep))
    sys.stderr.write("comp %d: %d valores, |V|>1e-6 hasta R=%.1f fm, "
                     "ahora %d validos + %d ceros\n"
                     % (comp, n, last*drf_old, len(toks[step-1::step][:npt]),
                        npt-len(toks[step-1::step][:npt])))
    out.append("%7d.%7d." % (npt, comp))
    for k in range(0, npt, 5):
        out.append(''.join('%16s' % t for t in keep[k:k+5]))
    i += 1 + nl
open(fout, 'w').write('\n'.join(out) + '\n')
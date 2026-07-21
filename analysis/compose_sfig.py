#!/usr/bin/env python3
# compose_sfig.py panelA.png panelB.png out.png
# Stack two per-locus splice panels vertically, add (a)/(b) labels + one-line descriptors.
import sys
from PIL import Image, ImageDraw, ImageFont
a_path, b_path, out = sys.argv[1], sys.argv[2], sys.argv[3]
A = Image.open(a_path).convert('RGB')
B = Image.open(b_path).convert('RGB')
W = max(A.width, B.width)
def pad(im):
    if im.width == W: return im
    c = Image.new('RGB', (W, im.height), 'white'); c.paste(im, ((W - im.width) // 2, 0)); return c
A, B = pad(A), pad(B)
LBL = 34   # top strip per panel for the (a)/(b) descriptor line
canvas = Image.new('RGB', (W, A.height + B.height + 2 * LBL), 'white')
canvas.paste(A, (0, LBL))
canvas.paste(B, (0, A.height + 2 * LBL))
d = ImageDraw.Draw(canvas)
def font(sz):
    for p in ['/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf',
              '/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf']:
        try: return ImageFont.truetype(p, sz)
        except Exception: pass
    return ImageFont.load_default()
fb = font(26); fs = font(19)
d.text((18, 6), '(a)', fill='black', font=fb)
d.text((70, 10), 'GM better-supported: the added exons are bridged by spliced RNA-seq reads (confirmed junctions)',
       fill='#1a1a1a', font=fs)
y2 = A.height + LBL + 6
d.text((18, y2), '(b)', fill='black', font=fb)
d.text((70, y2 + 4), 'GM worse-supported: the added ab initio exons carry no spliced reads or junctions (unconfirmed)',
       fill='#1a1a1a', font=fs)
canvas.save(out, dpi=(150, 150))
print('wrote', out, canvas.size)

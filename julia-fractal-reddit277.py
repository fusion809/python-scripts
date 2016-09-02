#!/usr/bin/env python3
# https://www.reddit.com/r/dailyprogrammer/comments/4v5h3u/20160729_challenge_277_hard_trippy_julia_fractals/
from PIL import Image

f = lambda q: q*q + complex(-0.221, -0.713)

w, h = 500, 400
img = Image.new('L', (w, h))
pix = img.load()

it_dep = 255
th = 2

n_min, n_max = -1, 1

for x, i in enumerate(x * ((n_max - n_min) / w) for x in range(int(n_min * (w/2)), int(n_max * (w/2)))):
    for y, r in enumerate(x * ((n_max - n_min) / h) for x in range(int(n_min * (h/2)), int(n_max * (h/2)))):
        g = 0
        z = complex(r, i)

        while abs(z) < th and g < it_dep:
            z = f(z)
            g += 1
        else:
            pix[x, y] = g

img.save('fractal.png')

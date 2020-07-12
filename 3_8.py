from fdtd3d import fdtd_gen
from math import sqrt, sin, pi
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

mi0 = 4e-7 * pi
eps0 = 8.854187817e-12

mi = 1.256627e-6
eps = 50e-12

sigma = 6e-8

delta = 0.01
delta_t = delta/3e8/sqrt(2)

i_max = j_max = 80

coefs = [[{
    'Ca': (1-sigma*delta/2/eps)/(1+sigma*delta/2/eps),
	'Cb': delta_t / (eps*delta)/(1+sigma*delta/2/eps),
	'Da': 1,
	'Db': delta_t / (mi*delta)
} for i in range(0, i_max)] for j in range(0, j_max)]

gen = fdtd_gen(i_max, delta, delta_t, coefs,
               f_meio=lambda t: sin(8e9*t))

Ez, Hx, Hy = next(gen)

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 8))
im = ax1.imshow(Ez, interpolation='none', aspect=1, vmin=-1, vmax=1)
im2 = ax2.imshow(Hx, interpolation='none', aspect=1, vmin=-0.002, vmax=0.002)
im3 = ax3.imshow(Hy, interpolation='none', aspect=1, vmin=-0.002, vmax=0.002)


def update(t):
	Ez, Hx, Hy = next(gen)
	im.set_array(Ez)
	im2.set_array(Hx)
	im3.set_array(Hy)
	return [im, im2, im3]


anim = animation.FuncAnimation(
	fig,
	update,
	frames=1240,
	interval=30
)
plt.show()
#anim.save('./onda.mp4')
#anim.save(os.path.dirname(os.path.realpath(__file__))+'\\onda.gif', writer='imagemagick')

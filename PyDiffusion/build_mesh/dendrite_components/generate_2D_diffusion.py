from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('TKAgg')
# matplotlib.use( "agg" )
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from math import *
from random import *

rand = np.random.uniform
randn = np.random.standard_normal
zeros = np.zeros
Mx = np.array # numpyArray = np.array([a,b,c,d])

#############################################################

NSteps = 100			# number of steps (loop time)
Ndots = 10				# number of particles
Scale = 1/10			# scale of model (life:model)
TimeStep = 1000			# time step (ms)
Delay = 0				# slow animation rate by X
DiffRateA = 0.5			# diffusion rate coefficient A
DiffRateB = 0.1			# diffusion rate coefficient B

Sc = Scale				# scale of model (life:model)
t = TimeStep/1000		# time step (ms)
dm = 2                  # dimensions
Da = DiffRateA*t/Sc		# Diffusion Rate A (D = L / 2d*t)
Db = DiffRateB*t/Sc		# Diffusion Rate B
Dr = Da/Db				# Ratio of Da:Ds (1/Ls)^2
Dn = Da/Dr				# new D after scaling L
k = sqrt(dm*Da)			# stdev of Ds step size distribution
L = sqrt(2*dm*Da)		# average diagonal (2D) step size
Lx = L/sqrt(2)          # average linear (1D) step size
Ls = 1/sqrt(Dr)			# scales Lx values for Dn
MSD = 2*dm*Da			# mean squared displacement

XYL = zeros((2,Ndots))		# XY particle locations
XYS = zeros((2,Ndots))		# XY step sizes
XYLp = zeros((2,NSteps))	# preallocate matrix for trace dot

# print XYL, XYS, XYLp

POLYSz = Mx([10, 15])		# size of polygon enclosure (XY in um)
POLYSz = POLYSz / Sc		# scale enclosures 
XWIDE = POLYSz[0] / 2		# half X enclosure size (will double below)
YHIGH = POLYSz[1] / 2		# half X enclosure size (will double below)
BOX = Mx([-1, 2, 2, 2]) / Sc 			# special area location [X Y W H]
BOARDER = Mx([-5, -7.5, 10, 15]) / Sc 	# [POLYSz(1)/2 POLYSz(2)/2 POLYSz(1) POLYSz(2)]

# print POLYSz, XWIDE, YHIGH, BOX


for Nt in range(0,NSteps):

	XYS = k * randn((2,Ndots))
	XYL = XYL+XYS

	XYLp[:,Nt] = XYL[:,0]

XLp = XYLp[0,:]
YLp = XYLp[1,:]

line, = plt.plot(XLp, YLp, '--', linewidth=2)
plt.show()

def update_line(num, data, line):
    line.set_data(data[...,:num])
    return line,

fig1 = plt.figure()
data = XYLp
l, = plt.plot([], [], 'r-')
plt.xlim(-sqrt(MSD*NSteps)*2, sqrt(MSD*NSteps)*2)
plt.ylim(-sqrt(MSD*NSteps)*2, sqrt(MSD*NSteps)*2)
plt.xlabel('x')
plt.title('test')
line_ani = animation.FuncAnimation(fig1, update_line, NSteps, fargs=(data, l),
    interval=50, blit=True)
# line_ani.save('lines.mp4')
plt.show()


# Uncomment to save images for animation
'''
for num in np.arange(NSteps):
	nvm = num+2
	plt.figure('abc' + str(num),figsize=(20,10),dpi=72)
	line, = plt.plot(XYLp[0,0:nvm], XYLp[1,0:nvm], '--', linewidth=2)
	plt.xlim(-sqrt(MSD*NSteps)*2, sqrt(MSD*NSteps)*2)
	plt.ylim(-sqrt(MSD*NSteps)*2, sqrt(MSD*NSteps)*2)
	plt.savefig("%d.png" % num,dpi=72)
	plt.close()
	print num
'''



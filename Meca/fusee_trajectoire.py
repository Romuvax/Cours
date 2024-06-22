#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


ForceFalcon = 7.607e6 # Newtons
m0 = 549e3 #mass in kg
alpha = 530e3/162. # mass loss in kg per second
#alpha = 0.

g = 9.81 # g
GravConstant = 6.67e-11 #N m^2/kg^2
RadiusEarth = 6400e3 #meters
MassEarth = 6e24 #kg

tmax = 162. #s
deltat = 0.1 # s

def grav_accel_of_altitude(height):
    """gravitational acceleration as a function of altitude in meters"""
    return GravConstant*MassEarth/(RadiusEarth+height)**2

def massnow(time):
    """"""
    return m0 - alpha*time


def StepVeloc(veloc, oldveloc, height, time):
    """"""
    newveloc  = oldveloc \
              + 2.*deltat \
              * (ForceFalcon \
              - massnow(time)*grav_accel_of_altitude(height) \
              + alpha * veloc)/massnow(time)
    return newveloc

def StepPosition(veloc, height, oldheight, time):
    """"""
    newheight = oldheight + 2*deltat * veloc
    return newheight

##############################################################33

nstep = int(tmax / deltat)

times = np.zeros(nstep)
heights = np.zeros(nstep)
velocs = np.zeros(nstep)

times[1] = deltat
velocs[1] = StepVeloc(velocs[0], velocs[0], heights[0], times[0])
heights[1] = StepPosition(velocs[0], heights[0], heights[0], times[0])

for istep in range(2,nstep):
    times[istep] = istep * deltat
    velocs[istep] = StepVeloc(velocs[istep-1], velocs[istep-2], heights[istep-1], times[istep-1])
    heights[istep] = StepPosition(velocs[istep-1], heights[istep-1], heights[istep-2], times[istep-1])
    
##############################################################33

#fig = plt.figure()
# create axis for present subplot
fig, (ax,ay) = plt.subplots(nrows=2,ncols=1,sharex=True)

#ax.set_xlim([0 , 2.*pi])
#ax.set_xticks(ax.get_xticks()[::1])
ax.set(ylabel='height(km)')
ax.plot(times,heights/1e3, label = "height aafo time")

#ay.set_xlim([0 , 2.*pi])
#ay.set_xticks(ax.get_xticks()[::1])
ay.set(xlabel='time(s)')
ay.set(ylabel='velocity(m/s)')
ay.plot(times,velocs)

# output to file
out_fig_name  =  "fusee_trajectoire.pdf"
fig.savefig( out_fig_name )
#plt.plot()

plt.close()

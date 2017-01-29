# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 18:09:45 2014

@author: Andreas

Edited on Wed Nov 26 
@Editor: David

Thanks: http://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
import reader as reader

f = open("NoReg_VarH_0.001.dat", 'r')

#zelfgeschreven klasse reader, neemt een geopende file als input.
#reader leest via update() de regels van de file een voor een uit en 
#slaat de relevante waarden op. Voor het uitlezen wordt gebruikgemaakt van 
#getters, hoewel denk ik ook rechtstreeks aan de variabelen geraakt kan worden.
rdr = reader.reader(f)
dt = 1./30 # 30 fps

energyP = []
timeP = []


# set up figure and animation
gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
fig = plt.figure(facecolor="beige")


ax = plt.subplot(gs[0])
ax.set_ylim([-2,2])
ax.set_xlim([-2,2])
ax.set_xlabel('x ')
ax.set_ylabel('y')


axE = plt.subplot(gs[1])
axE.set_yscale('log')
axE.set_ylim([1e-16,0.1])
axE.set_xlim([0,900])
axE.set_xlabel('Time')
axE.set_ylabel('Fractional Energy Error')


plt.tight_layout()


# particles hold the locations of the particles

particle1, = ax.plot([], [], color='darkorange',marker='o', ms=6)
particle2, = ax.plot([], [], color='cyan',marker='o', ms=6)
particle3, = ax.plot([], [], color='limegreen',marker='o', ms=6)
energy, = axE.plot([], [], color='darkcyan', linewidth=2, ms=6)

#strings die de waarden tonen

time_text = ax.text(0.05, 0.85, '', transform=ax.transAxes)

def init():
    """initialize animation"""
    particle1.set_data([], [])
    particle2.set_data([], [])
    particle3.set_data([], [])
    energy.set_data([], [])
    time_text.set_text('')
    return particle1, particle2, particle3, energy, time_text

def animate(i):
    """perform animation step"""
    global rdr, dt, energyP, timeP
    rdr.update()

    energyP.append(rdr.getEnergy())
    timeP.append(rdr.getT())
    dataP = [timeP, energyP]

    particle1.set_data(rdr.getPositions(0))
    particle2.set_data(rdr.getPositions(1))
    particle3.set_data(rdr.getPositions(2))
    energy.set_data(dataP)
    
    time_text.set_text('time = %.1f' % rdr.getT())
    return particle1, particle2, particle3, energy, time_text


# choose the interval based on dt and the time to animate one step.
from time import time
t0 = time()
animate(0)
t1 = time()
#zet op 1000 voor goed tempo! 
interval = 100 * dt - (t1 - t0)

#animatie wordt gemaakt. animate is de functie die blijft aangeroepen worden 
#tot oneindig. Geen manier gevonden om te stoppen. Als de file volledig uitgelezen
#is, geeft update voor alle waarden 0 terug. (zie reader). 'frames' is belangrijk
#bij het opslaan van de animatie: de animatie wordt maar opgeslagen voor de 
#hoeveelheid frames die hier vermeld staan!

ani = animation.FuncAnimation(fig, animate, interval=interval, frames=500, blit=True, init_func=init)


# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html

#schrijven van nieuwe ffmpegwriter, voor fout in broncode

mywriter = animation.FFMpegWriter(30, bitrate=-1)

#uncomment volgende lijn voor het opslaan van de file.

ani.save('BurrauNoReg.mp4', writer=mywriter, fps=30, extra_args=['-vcodec', 'libx264'])

#sluiten van f.

f.close

plt.show()

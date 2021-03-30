from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.animation as animation
import numpy as np
import os
from main import ttab, ttab1
from Viriel import R_final, tvir1

os.chdir("C:/Users/Lucas/PycharmProjects/Projet_Num")

fig = plt.figure(figsize=(15, 15))
ax = fig.gca(projection='3d')


# Fonction qui va réactualiser le rayon R au temps t_n
def update(num, x, y, z, plot):
    plot[0].remove()
    r = R_final[num]
    x_ = x * r
    y_ = y * r
    z_ = z * r
    title = 'Rayon au temps t = ' + str(round(ttab1[num], 2)) + "  milliards d'années"
    ax.set_title(title, fontsize=15)
    plot[0] = ax.plot_surface(x_, y_, z_, cmap="hot")
    if ttab[num] > tvir1:
        fake2Dline = matplotlib.lines.Line2D([0], [0], linestyle="none", c='orange', marker='o')
        plot[0] = ax.plot_surface(x_, y_, z_, cmap="hot")
        ax.legend([fake2Dline], ["Rayon de virialisation atteint"], fontsize=12, numpoints=1)


# Coordonnées de la sphère

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

x = 2 * np.outer(np.cos(u), np.sin(v))
y = 2 * np.outer(np.sin(u), np.sin(v))
z = 2 * np.outer(np.ones(np.size(u)), np.cos(v))

# Délimitation des axes
ax.set_xlim(- R_final[-1] * 5, R_final[-1] * 5)
ax.set_ylim(- R_final[-1] * 5, R_final[-1] * 5)
ax.set_zlim(- R_final[-1] * 5, R_final[-1] * 5)
ax.set_axis_off()
plot = [ax.plot_surface(x, y, z)]
ax.set_facecolor('k')

matplotlib.rcParams['animation.ffmpeg_path'] = r'C:/Users/Lucas/ffmpeg/bin/ffmpeg.exe'

anim = animation.FuncAnimation(fig, update, 450, fargs=(x, y, z, plot), interval=1e-10)
# POUR REALISER L'ANIMATION,REPASSER à MAXIMUM 5000 INTERVALLES DE TEMPS (N=5000)
# ET 600 FRAMES,SINON L'OPERATION DEVIENT TROP LOURDE


writervideo = animation.FFMpegWriter(fps=35)
# anim.save("animation_finale.mp4", writer=writervideo)

plt.show()

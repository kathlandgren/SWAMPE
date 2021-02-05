from matplotlib import animation
from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap


fig, ax = plt.subplots()

# set up map projection
m = Basemap(projection='nsper',lon_0=-0,lat_0=90)
m.drawcoastlines()
m.drawparallels(np.arange(0.,180.,30.))
m.drawmeridians(np.arange(0.,360.,60.))

# some 2D geo arrays to plot (time,lat,lon)
data = np.random.random_sample((20,90,360))
lat = np.arange(len(data[0,:,0]))
lon = np.arange(len(data[0,0,:]))
lons,lats = np.meshgrid(lon,lat)

# ims is a list of lists, each row is a list of artists to draw in the
# current frame; here we are animating three artists, the contour and 2 
# annotatons (title), in each frame
ims = []
for i in range(len(data[:,0,0])):
    im = m.contourf(lons,lats,data[i,:,:],latlon=True)
    add_arts = im.collections
    text = 'title={0!r}'.format(i)
    te = ax.text(90, 90, text)
    an = ax.annotate(text, xy=(0.45, 1.05), xycoords='axes fraction')
    ims.append(add_arts + [te,an])

ani = animation.ArtistAnimation(fig, ims)
## If you have ffmpeg you can save the animation by uncommenting 
## the following 2 lines
# FFwriter = animation.FFMpegWriter()
# ani.save('basic_animation.mp4', writer = FFwriter)
plt.show()

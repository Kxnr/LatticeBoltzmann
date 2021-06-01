import numpy as np
import matplotlib.pyplot as mp
from matplotlib import cm
import sys
import os


def curl(image):
    # naiive method is super slow!
    # TODO: maybe try to optimize with convolution
    c, w, h = image.shape
    vorticity = np.zeros((w, h))
    for x in range(w):
        for y in range(h):
            vorticity[x, y] = (image[1, (x+1) % w, y] -
                               image[1, (x-1) % w, y]) - \
                (image[0, x, (y+1) % h] -
                 image[0, x, (y-1) % h])

    vorticity /= 2.
    return vorticity


if __name__ == "__main__":
    prefix = sys.argv[1]
    num_frames = int(sys.argv[2])

    vmax = None
    vmin = None
    cmax = None
    cmin = None

    for i in range(num_frames):
        f = np.genfromtxt(prefix+repr(i)+".txt", skip_header=2)
        f = f.reshape((2, -1, f.shape[-1]))
        u = np.sum(f**2, axis=0)
        c = curl(f)

        # get color bar limits
        if vmax is None and vmin is None:
            vmax = np.amax(u)
            vmin = np.amin(u)

            vmax += .1*np.abs(vmax)
            vmin -= .1*np.abs(vmin)

        if cmax is None and cmin is None:
            cmax = np.amax(c)
            cmin = np.amin(c)

            cmax += .1*np.abs(cmax)
            cmin -= .1*np.abs(cmin)

        mp.clf()
        
        mp.imshow(u.transpose(), cmap=cm.plasma, vmin=vmin, vmax=vmax)
        mp.colorbar()
        mp.savefig("velocity."+str(i)+".png")

        mp.clf()
        mp.imshow(u.transpose(), cmap=cm.plasma)
        mp.colorbar()
        mp.savefig("velocity_nonorm."+str(i)+".png")

        mp.clf()
        mp.imshow(c.transpose(), cmap=cm.seismic, vmin=cmin, vmax=cmax)
        mp.colorbar()
        mp.savefig("curl."+str(i)+".png")

        mp.clf()
        mp.imshow(c.transpose(), cmap=cm.seismic)
        mp.colorbar()
        mp.savefig("curl_nonorm."+str(i)+".png")

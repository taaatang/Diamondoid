import numpy as np
from numba import njit

@njit
def projectTo2D(coords):
    L = coords.shape[0]
    coords2d = np.zeros((L, 2))
    zmin = np.min(coords[:,2])
    zmax = np.max(coords[:,2])
    zd = 2.0 / (zmax - zmin)
    for i in range(L):
        x = coords[i, 0]
        y = coords[i, 1]
        z = coords[i, 2]
        xsign = 1
        ysign = 1
        d = z % 4
        if d>=2:
            xsign = -1
        if d==1 or d==3:
            ysign = -1
        offset = zd * (z - zmin) - 1.0
        coords2d[i, 0] = x + xsign * offset
        coords2d[i, 1] = y + ysign * offset
    return coords2d

@njit
def findNN(i, coords):
    indices = np.ones(4, dtype=np.int32) * (-1)
    count = 0
    rmin = 0.9*np.sqrt(3)
    rmax = 1.1*np.sqrt(3)
    for j in range(i + 1, coords.shape[0]):
        d = np.linalg.norm(1.0*(coords[i] - coords[j]))
        if d>rmin and d<rmax:
            indices[count] = j
            count += 1
            if count >= 4:
                break
    return indices
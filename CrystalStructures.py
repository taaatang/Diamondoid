import numpy as np
from mayavi.mlab import *





# One unit-cell of diamond contains 8 atoms,
# coordates in unit of 3.566700 Angstrom (the lattice constant of diamond)
# C-C bond length is 1.54443 Angstrom

DiamondCoords = np.zeros((8,3))
DiamondCoords[0,:] = [0.000000,   0.000000,   0.000000]
DiamondCoords[1,:] = [0.000000,   0.500000,   0.500000]
DiamondCoords[2,:] = [0.500000,   0.500000,   0.000000]
DiamondCoords[3,:] = [0.500000,   0.000000,   0.500000]
DiamondCoords[4,:] = [0.750000,   0.250000,   0.750000]
DiamondCoords[5,:] = [0.250000,   0.250000,   0.250000]
DiamondCoords[6,:] = [0.250000,   0.750000,   0.750000]
DiamondCoords[7,:] = [0.750000,   0.750000,   0.250000]


nAdamantaneCoords = 10
AdamantaneCoords = np.zeros((nAdamantaneCoords,3))
AdamantaneCoords[0,:] = [0.000000,   0.500000,   0.500000]
AdamantaneCoords[1,:] = [1.000000,   0.500000,   0.500000]
AdamantaneCoords[2,:] = [0.500000,   0.500000,   0.000000]
AdamantaneCoords[3,:] = [0.500000,   0.500000,   1.000000]
AdamantaneCoords[4,:] = [0.500000,   0.000000,   0.500000]
AdamantaneCoords[5,:] = [0.500000,   1.000000,   0.500000]
AdamantaneCoords[6,:] = [0.750000,   0.250000,   0.750000]
AdamantaneCoords[7,:] = [0.250000,   0.250000,   0.250000]
AdamantaneCoords[8,:] = [0.250000,   0.750000,   0.750000]
AdamantaneCoords[9,:] = [0.750000,   0.750000,   0.250000]

# {0.75, 0.25, -0.25},\
# {0.25, 0.75, -0.25},\
# {0.0, 0.0, 0.0},\
# {1.0, 1.0, 0.0},\
# {0.75, -0.25, 0.25},\
# {0.25, -0.25, 0.75},\
# {1.25, 0.25, 0.25},\
# {1.25, 0.75, 0.75},\
# {-0.25, 0.75, 0.25},\
# {-0.25, 0.25, 0.75},\
# {0.25, 1.25, 0.25},\
# {0.75, 1.25, 0.75},\
# {1.0, 0.0, 1.0},\
# {0.0, 1.0, 1.0},\
# {0.25, 0.25, 1.25},\
# {0.75, 0.75, 1.25}


# potential:

nTriamantaneCoords = 18
TriamantaneCoords = np.zeros((nTriamantaneCoords,3))
TriamantaneCoords[0,:] = [-1.000000,   0.000000,   1.000000]
TriamantaneCoords[1,:] = [0.000000,   0.000000,   0.000000]
TriamantaneCoords[2,:] = [-1.000000,   0.500000,   0.500000]
TriamantaneCoords[3,:] = [0.000000,   0.500000,   0.500000]
TriamantaneCoords[4,:] = [-0.500000,   0.500000,   0.000000]
TriamantaneCoords[5,:] = [-0.500000,   0.500000,   1.000000]
TriamantaneCoords[6,:] = [-0.500000,   0.000000,   0.500000]
TriamantaneCoords[7,:] = [-0.500000,   1.000000,   0.500000]
TriamantaneCoords[8,:] = [-0.250000,   0.250000,   0.750000]
TriamantaneCoords[9,:] = [-0.750000,   0.250000,   0.250000]
TriamantaneCoords[10,:] = [0.250000,   0.250000,   0.250000]
TriamantaneCoords[11,:] = [-0.750000,   0.750000,   0.750000]
TriamantaneCoords[12,:] = [-0.250000,   0.750000,   0.250000]
TriamantaneCoords[13,:] = [-1.250000,   0.250000,   0.750000]
TriamantaneCoords[14,:] = [-0.750000,  -0.250000,   0.750000]
TriamantaneCoords[15,:] = [-0.250000,  -0.250000,   0.250000]
TriamantaneCoords[16,:] = [-0.250000,   0.250000,  -0.250000]
TriamantaneCoords[17,:] = [-0.75, 0.25, 1.25]


nTetramantane121Coords = 22
Tetramantane121Coords = np.zeros((nTetramantane121Coords,3))
Tetramantane121Coords[0,:] = [-1.000000,   0.000000,   1.000000]
Tetramantane121Coords[1,:] = [0.000000,   0.000000,   0.000000]
Tetramantane121Coords[2,:] = [0.000000,   1.000000,   0.000000]
Tetramantane121Coords[3,:] = [-1.000000,   0.500000,   0.500000]
Tetramantane121Coords[4,:] = [0.000000,   0.500000,   0.500000]
Tetramantane121Coords[5,:] = [-0.500000,   0.500000,   0.000000]
Tetramantane121Coords[6,:] = [-0.500000,   0.500000,   1.000000]
Tetramantane121Coords[7,:] = [0.500000,   0.500000,   0.000000]
Tetramantane121Coords[8,:] = [-0.500000,   0.000000,   0.500000]
Tetramantane121Coords[9,:] = [-0.500000,   1.000000,   0.500000]
Tetramantane121Coords[10,:] = [-0.250000,   0.250000,   0.750000]
Tetramantane121Coords[11,:] = [-0.750000,   0.250000,   0.250000]
Tetramantane121Coords[12,:] = [0.250000,   0.250000,   0.250000]
Tetramantane121Coords[13,:] = [-0.750000,   0.750000,   0.750000]
Tetramantane121Coords[14,:] = [-0.250000,   0.750000,   0.250000]
Tetramantane121Coords[15,:] = [-1.250000,   0.250000,   0.750000]
Tetramantane121Coords[16,:] = [-0.750000,  -0.250000,   0.750000]
Tetramantane121Coords[17,:] = [-0.250000,  -0.250000,   0.250000]
Tetramantane121Coords[18,:] = [-0.250000,   0.250000,  -0.250000]
Tetramantane121Coords[19,:] = [0.250000,   0.750000,  -0.250000]
Tetramantane121Coords[20,:] = [0.0,  0.5,  -0.5]
Tetramantane121Coords[21,:] = [-0.75, 0.25, 1.25]


nPentamantane1234Coords = 26
Pentamantane1234Coords = np.zeros((nPentamantane1234Coords,3))
Pentamantane1234Coords[0,:] = [1.  ,  1.  ,  1.  ]
Pentamantane1234Coords[1,:] = [ 0.  ,  1.  ,  1.  ]
Pentamantane1234Coords[2,:] = [ 1.  ,  0.  ,  1.  ]
Pentamantane1234Coords[3,:] = [ 1.  ,  1.  ,  0.  ]
Pentamantane1234Coords[4,:] = [ 1.  ,  0.5 ,  0.5 ]
Pentamantane1234Coords[5,:] = [ 1.  ,  0.5 ,  1.5 ]
Pentamantane1234Coords[6,:] = [ 1.  ,  1.5 ,  0.5 ]
Pentamantane1234Coords[7,:] = [ 0.5 ,  0.5 ,  1.  ]
Pentamantane1234Coords[8,:] = [ 0.5 ,  1.5 ,  1.  ]
Pentamantane1234Coords[9,:] = [ 1.5 ,  0.5 ,  1.  ]
Pentamantane1234Coords[10,:] = [ 0.5 ,  1.  ,  0.5 ]
Pentamantane1234Coords[11,:] = [ 0.5 ,  1.  ,  1.5 ]
Pentamantane1234Coords[12,:] = [ 1.5 ,  1.  ,  0.5 ]
Pentamantane1234Coords[13,:] = [ 0.75,  0.25,  0.75]
Pentamantane1234Coords[14,:] = [ 0.75,  1.25,  0.75]
Pentamantane1234Coords[15,:] = [ 0.25,  1.25,  1.25]
Pentamantane1234Coords[16,:] = [ 1.25,  0.25,  1.25]
Pentamantane1234Coords[17,:] = [ 1.25,  1.25,  0.25]
Pentamantane1234Coords[18,:] = [ 0.25,  0.75,  0.75]
Pentamantane1234Coords[19,:] = [ 1.25,  0.75,  0.75]
Pentamantane1234Coords[20,:] = [ 0.75,  0.75,  0.25]
Pentamantane1234Coords[21,:] = [ 0.75,  0.75,  1.25]
Pentamantane1234Coords[22,:] = [ 0.25,  0.25,  0.25]
Pentamantane1234Coords[23,:] = [ 0.5 ,  0.  ,  0.5 ]
Pentamantane1234Coords[24,:] = [ 0.  ,  0.5 ,  0.5 ]
Pentamantane1234Coords[25,:] = [ 0.5 ,  0.5 ,  0.  ]


# Visualization


for ix in range(3):
    for iy in range(3):
        for iz in range(3):
            points3d(AdamantaneCoords[:,0] + ix, AdamantaneCoords[:,1] + iy, \
                     AdamantaneCoords[:,2] + iz, resolution = 16, color = (1,1,1), \
                     scale_factor=0.2)

#points3d(Tetramantane121Coords[:,0], Tetramantane121Coords[:,1], \
#         Tetramantane121Coords[:,2], resolution = 16, scale_factor=0.2)

#points3d(TriamantaneCoords[:,0], TriamantaneCoords[:,1], \
#         TriamantaneCoords[:,2], resolution = 16, scale_factor=0.2)

#points3d(AdamantaneCoords[:,0], AdamantaneCoords[:,1], \
#         AdamantaneCoords[:,2], resolution = 16, color = (0,1,1), scale_factor=0.2)

points3d(Pentamantane1234Coords[:,0], Pentamantane1234Coords[:,1], \
         Pentamantane1234Coords[:,2], resolution = 16, color = (0,1,1), scale_factor=0.2)

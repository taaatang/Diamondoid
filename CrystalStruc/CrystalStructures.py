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
AdamantaneCoords = np.zeros((10,3))
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



# Visualization

points3d(Tetramantane121Coords[:,0], Tetramantane121Coords[:,1], \
         Tetramantane121Coords[:,2], resolution = 16, scale_factor=0.2)
         
         
         
         
         
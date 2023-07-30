# Use VMD to visualize molecular structures
In each data directory, data files are organized into sub directories: moleculeName/carbonAtomNumber/, each contains 10 data files from 10 independent simulations.

To visualize one data file:
1. open VMD, choose “File" -> "New Molecule"
2. in the “Filename" filed, click “Browse", then navigate to the data file and then click “Open"
3. in “Determine file type” filed, choose VASP_POSCAR, then click “Load”

Now you can play around with the loaded structure. You can go to “Graphics”->”Representations”, and change the “Drawing Method” to “CPK” to visualize C atoms and C-C bonds.
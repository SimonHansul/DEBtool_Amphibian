# AmP fitting code for amphibian metamorphosis model

**Author:** C. Romoli - Ecological modelling dept. - ibacon GmbH

**Date:** 09.08.2024

# Description of the folder

The folder contains all the necessary files to fit the amphibian metamorphosis module integrating with the "Add my Pet" (AmP) routines (more info on AmP [here](https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/index.html)).

In the subfolder `MATLAB_library`, there are the internal functions needed to solve the tadpole development during the metamorphosis process. The files in this folders are called by the `predict_<species>.m` files contained in the code folders specific for every species.

For some specific datasets relative to *Dryohpytes versicolor* and *Xenopus laevis*, particular versions of the general files contained in the `MATLAB_library` needed to be used. These dedicated versions are stored in the relative folder of the species.

To run the AmP fitting, the [DEBtool_M](https://github.com/add-my-pet/DEBtool_M) and [AmPtool](https://github.com/add-my-pet/AmPtool) are needed.

Because *Dryohpytes versicolor* is a species that is currently not present in the AmP database, it needs to be added. This can be achieved by copying the files `Hylidae.txt` and `Dryophytes.txt` in the `AmPtool/taxa` folder.

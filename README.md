# RDLA
Ramalina siliquosa diffusion limited aggregation

RDLA also requires the Eigen 
Eigen/Dense

and boost 
boost/function_output_iterator
boost/geometry/geometry

Installation:

Win10
Enable the Windows Subsystem for Linux (WSL) and install a distro: https://docs.microsoft.com/en-us/windows/wsl/install-win10
Use the Linux bash to run the install commands for Unix below.

Unix/Linux
Install make, boost and g++

sudo apt install make libboost-dev g++

Compile the project using the Makefile

make

Run the compiled application

./dlaf

# RDLA sheathing, sampling, and denoising

sheath-sample-denoise folder contains Python files to give the RDLA simulation volume, simulate laser scanning on the RDLA, and then denoise the lichen simulation to reveal the topography of the mesh below the simuation. This requires the pyCloudCompare Python Wrapper to run the CloudCompare CLI (included), pyMeshLab, and open3D libraries. Alternatively, download jupyter and run it in jupyter notebook (https://jupyter.org/), which includes these libraries 

Acknowledgements:
Jason Webb for the base DLA implementation: https://web.archive.org/web/20210822065831/https://github.com/jasonwebb/dlaf-with-models

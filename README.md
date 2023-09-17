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

Acknowledgements:
Jason Webb for the base DLA implementation: https://web.archive.org/web/20210822065831/https://github.com/jasonwebb/dlaf-with-models

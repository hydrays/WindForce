# WindForce
Powered by bigdata, this code is for estimation, prediction and visualization of air pollution.

# Pre-request
GDAL
libDAI

# Install
```
mkdir build
cd build
cmake ..
make
```

# Test with synthetic data
Consider a Ising model with Mdim x Ndim grids. We want to estimate the model parameters from partial observations.
To test the algorithm, three programs, namely the GB, BP and II, are used.

GB is used to generate synthetic data, the parameters are set in the code and will be stored in "parameters.txt". 
First, GP will generate the full state vector using parameters provides. Then we may want to cover some of the grids (missing data) and test the algorithm. The generated data will be stored in the "data" folder.

BP is used to fill the missing data using parameters listed in "result.txt". 
Note that initially these parameters are guessed. But after many interations of BP+II, we hope to learn the true parameters used to generate the data. The recovered data will be stored in the "evidence" folder.
For unknown grids, we approximate their s0 using nearest neighbour.

II is used to solve the inverse problem.
Note that II only estimate s0 on the known site. The unknown site will be printed as 0, but actually it is given by the approximation we just described.
Currently we don't need to change this file.

Notes:
1. In both "data" and "evidence" folder, make sure the "filelist.txt" contains all the *.tif files.
2. There are two places need to be modified in generate_data.cpp and one plance to be modified in belief_propagation.cpp. These places are marked as "modify".

# Contact:
Yucheng Hu
Zhou Pei-yuan Center for Applied Mathematics, Tsinghua University, Beijing, China, 100084.
huyc@tsinghua.edu.cn

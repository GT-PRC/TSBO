# Two-Stage Bayesian Optimization (TSBO)
Copyright (c) 2020 <br />
3D Systems Packaging Research Center (PRC) <br />
Georgia Insitute of Technology <br />

This material is based on work supported by NSF I/UCRC Center for Advanced Electronics Through Machine Learning (CAEML).
For questions and queries, please contact: htorun3@gatech.edu

Please cite our paper if you use the code: <br />
H. M. Torun, M. Swaminathan, A. Kavungal Davis and M. L. F. Bellaredj, <br />
["A Global Bayesian Optimization Algorithm and Its Application to Integrated System Design"](https://ieeexplore.ieee.org/document/8253829) <br />
in IEEE Transactions on Very Large Scale Integration (VLSI) Systems, vol. 26, no. 4, pp. 792-802, April 2018.
## Scripts
###### TSBO.m <br />
The main code including implementation of TSBO, definition of objective objective function and parameters of TSBO.
The code is commented to explain implementation details and choice of default hyperparameters. Additional details 
on choice of parameters can be found in TVLSI paper Section III.D

###### splitregion.m <br />
Used in 1st stage of TSBO for hierarchical partioning tree. The script takes a D-dimensional sample space as input in the form of Dx2 matrix, 
where columns represent min and max of sample space for each parameter and divides it into 2^D new regions. The output is Dx2x(2^D) matrix,
where each Dx2 matrix contains a partition of the inputed Dx2 sample space.

###### splitregion33.m <br />
Used in 2nd stage of TSBO. Same as the "splitregion.m", but divides the input sample space into 3 new regions along the largest vector.
The output is Dx2x3 matrix. 


###### "test_function" <br />
This folder contains sample objective functions that can be used to test the algorithm. There are multiple functions defined in the folder. These are "Optimization Challenge Functions", each function has different response surface and many local optima.

## System Requirements
The code is tested using Matlab R2018b and R2019a, and has dependency to "Statistics and Machine Learning Toolbox"

Copyright (c) 2018 Hakki Mert Torun, Georgia Insitute of Technology.

Two-Stage Bayesian Optimization (TSBO) for maximizing a black-box function.

This material is based on work supported by NSF I/UCRC Center for Advanced Electronics Through Machine Learning (CAEML).
For questions and queries, please contact: htorun3@gatech.edu, 404-528-5943

Please cite our paper if you use the code:
H. M. Torun, M. Swaminathan, A. Kavungal Davis and M. L. F. Bellaredj, 
"A Global Bayesian Optimization Algorithm and Its Application to Integrated System Design," 
in IEEE Transactions on Very Large Scale Integration (VLSI) Systems, vol. 26, no. 4, pp. 792-802, April 2018.
---------------------------------------------------------------------------------------------------------------------------
We tested the code on Matlab 2017a, 2017b, 2018a (Windows 7 & 10 versions).
The folder should contain following scripts:

--TSBO.m
The main code including implementation of TSBO, definition of objective objective function and parameters of TSBO.
The code is commented to explain implementation details and choice of default hyperparameters. Additional details 
on choice of parameters can be found in TVLSI paper Section III.D

--splitregion.m
Used in 1st stage of TSBO for hierarchical partioning tree. The script takes a D-dimensional sample space as input in the form of Dx2 matrix, 
where columns represent min and max of sample space for each parameter and divides it into 2^D new regions. The output is Dx2x(2^D) matrix,
where each Dx2 matrix contains a partition of the inputed Dx2 sample space.

--splitregion33.m
Used in 2nd stage of TSBO. Same as the "splitregion.m", but divides the input sample space into 3 new regions along the largest vector.
The output is Dx2x3 matrix. 


--Folder: "test_function"
This contains sample objective functions that can be used to test the algorithm. There are multiple functions defined in the folder.
These are "Optimization Challenge Functions", each function has different response surface and many local optima.
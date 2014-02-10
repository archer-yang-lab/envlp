% AUXILIARY
% MATLAB Version 8.0 (R2012b) 12-Feb-2014
%
% Files
%   center           - Subtract the mean of each column
%   Contr            - Compute the contraction matrix of dimension r
%   Expan            - Compute the expansion matrix of dimension r
%   fit_OLS          - Fit multivariate linear regression
%   get_combeig      - This function gives an array (A) of nchoosek(p,u) matrices each one being a permutation of the columns of matrix V
%   get_envelope     - Construct the envelope subspace using a sequential algorithm
%   get_Init         - Generate starting value for the envelope subspace
%   get_Init4envmean - Generate starting value for the envelope subspace in estimating the multivariate mean
%   get_Init4henv    - Generate starting value for the heteroscedastic envelope subspace
%   grams            - Gram-Schmidt orthogonalization of the columns of A
%   Kpd              - Compute the communication matrix Kpd
%   Lmatrix          - Extract the 2nd to the last diagonal element of a matrix into a vector
%   make_opts        - Make optional input parameters for running the sg_min package
%   make_parameter   - Compute summary statistics from the data
%   MBoxtest         - Multivariate Statistical Testing for the Homogeneity of Covariance Matrices by the Box's M
%   mtest            - Perform Box's M test to check the homogeneity of the covariance matrices
%   nulbasis         - Compute nulbasis basis for nullspace

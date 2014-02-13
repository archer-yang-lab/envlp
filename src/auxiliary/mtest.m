%% mtest
% Perform Box's M test to check the homogeneity of the covariance matrices. 

%% Syntax
%         TestOutput = mtest(X, Y, alpha)
%
%% Input
%
% *X*: Group indicators. A matrix with n rows.  X can only have p unique
%  rows, where p is the number of groups. For example, if there 
% are two groups, X can only have 2 different kinds of rows, such as (0, 1)
% and (1, 0), or (1, 0, 10) and (0, 5, 6).  The number of columns is not
% restricted, as long as X only has p unique rows.
%
% *Y*: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
%
% *alpha*: Significance level for testing.  A real number between 0 and 1,
% often taken at 0.05 or 0.01.

%% Output
% 
% *TestOutput*: A list containing the Box's M statistic, the approximation
% test statistic, degrees of freedom for the approximation statistic test,
% and the p-value.  At the same time, a table is printed out.
% 
% * TestOutput.mStatistic: The Box's M statistic. A real number. 
% * TestOutput.approxStatistic: The approximation test statistic.
% * TestOutput.df: The degrees of freedom of the approximation statistic
% test.  A positive integer. 
% * TestOutput.pValue: p-value of the test.  A real number in [0, 1].

%% Description
% This function performs the Box's M test for homegeneity of the covariance
% matrices for different groups, indicated by X. If the groups sample-size 
% is at least 20 (sufficiently large), Box's M test takes a Chi-square 
% approximation; otherwise it takes an F approximation.

%% References
% The codes are implemented based on 
% Trujillo-Ortiz, A., R. Hernandez-Walls, K. Castro-Morales, A. Espinoza-Tenorio, A. Guia-Ramirez
% and R. Carmona-Pina. (2002). MBoxtest: Multivariate Statistical Testing for the Homogeneity of 
% Covariance Matrices by the Box's M. A MATLAB file. [WWW document]. URL http://www.mathworks.com/
% matlabcentral/fileexchange/loadFile.do?objectId=2733&objectType=FILE

%% Example
%         load waterstrider.mat
%         alpha = 0.01;
%         TestOutput = mtest(X, Y, alpha);

function TestOutput = mtest(X, Y, alpha)

[n, r] = size(Y);
[Xs, ~] = sortrows(X);
xMarker = ones(n, 1);

temp = Xs(1, :);
p = 1;
for i = 1 : n

    if prod(double(Xs(i, :) == temp)) == 0            
        temp = Xs(i, :);
        p=p+1;
    end
    xMarker(i) = p;
    
end 

Z = zeros(n, r + 1);
Z(:, 1) = xMarker;
Z(:, 2 : end) = Y;
[MB, X2, v, P] = MBoxtest(Z, alpha);
TestOutput.mStatistic = MB;
TestOutput.approxStatistic = X2;
TestOutput.df = v;
TestOutput.pValue = P;

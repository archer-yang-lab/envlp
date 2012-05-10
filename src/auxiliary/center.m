%% center
% Substract the mean of each column

%% Syntax
%         XC = center(X)
%
%% Input
%
% *X*: A matrix or a column vector. 
%
%% Output
%
% *XC*: A matrix or a column vector with the mean for each column equal to
% 0.

%% Description
% This function centerizes a matrix or a vector, by substracting each
% column by its column mean.

%% Example
%
%         a = [1 2 3; 4 5 6];
% 
%         center(a)

function XC = center(X);

n = size(X, 1);

XC = X - ones(n, 1) * mean(X);


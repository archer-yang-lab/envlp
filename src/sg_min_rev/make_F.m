%% make_F
% Generic function to generate the objective function F.
 
%% Syntax
% F = make_F(fun_method_handle, FParameters)
% 
% Input
%
% * fun_method_handle: A specific model objective function.
% * FParameters: A structure that contains data parameters as input for
% the function fun_method_handle.
%
% Output
%
% * F: The generic objective function for computing the inner envelope subspace.
%

%% Description
%
% Generic function to generate the objective function F. The function first 
% sets a handle to the specific model function and fixes the data parameters from 
% the sample needed for its computation. The handle fixed with those 
% parameters is then evaluated at a given value for argument W. A generic objective
% function F is returned.

function F = make_F(fun_method_handle, FParameters)

F = @tmp;
    function fval = tmp(W, varargin)
        fval = fun_method_handle(W, FParameters, varargin{:});
    end
end


%% make_dF
% Generic function to generate the derivative function of the objective function F.
 
%% Syntax
% dF = make_dF(dfun_method_handle, FParameters)
% 
%% Input
%
% * dfun_method_handle: A specific model derivative function of the objective function.
% * FParameters: A structure that contains data parameters as input for
% the function dfun_method_handle.
%
%% Output
%
% * dF: The generic derivative function of the objective function for computing the envelope subspace.
%

%% Description
%
% Generic function to generate the derivative function of the objective function F. 
% The function first sets a handle to the specific model function and fixes the data 
% parameters from the sample needed for its computation. The handle fixed with those 
% parameters is then evaluated at a given value for argument W. A generic derivative 
% function dF is returned.

function dF = make_dF(dfun_method_handle, FParameters)

dF = @tmp;
    function dfval = tmp(W)
        dfval = dfun_method_handle(W, FParameters);
    end
end


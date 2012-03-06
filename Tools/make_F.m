function f = make_F(fun_method_handle,FParameters)
%
% f = F(fun_handle,FParameters)
%
% Generic function to compute the objective function F. The function first 
% sets a handle to the specific model function and fixes the parameters from 
% the sample needed for its computation. The handle fixed with those 
% parameters is then evaluated at a given value for argument W.
%
%   fun_handle: handle to the objective function for a specific model.
%   FParameters: fixed parameters needed to evaluate the objective function.
%   W: argument for computing the derivative.
% =========================================================================

f = @tmp;
    function fval = tmp(W,varargin)
        fval = fun_method_handle(W,FParameters,varargin{:});
    end
end


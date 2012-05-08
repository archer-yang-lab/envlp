function dF = make_dF(dfun_method_handle,FParameters)
%
% diff = dF(dfun_handle,FParameters)
%
% Generic function to compute the derivative of F. The function first sets
% a handle to the specific derivative and fixes the parameters from the
% sample needed for its computation. The handle fixed with those parameters
% is then evaluated at a given value for argument W.
%
%   dfun_handle: handle to the function derivative for a specific model.
%   FParameters: fixed parameters needed to compute the derivative.
%   W: argument for computing the derivative.
% =========================================================================
dF = @tmp;
    function dfval = tmp(W)
        dfval = dfun_method_handle(W,FParameters);
    end
end


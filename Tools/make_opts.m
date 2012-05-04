%% make_opts
% Make optional input parameters for running the sg_min package.

%% Usage
% opts = make_opts(opts)
%
% Input 
%
% A list containing optional input parameter for sg_min.m specified by users.  One or several (even all)
% fields could be empty.
%
% * opts.maxIter: Maximum number of iterations.  Default value: 300.
% * opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * opts.verbose: Flag for print out output, logical 0 or 1. Default value: 0.
%
% Output:
% 
% A list containing optional input parameter for sg_min.m, specified by
% users or the default values are used.
%
% * opts.maxIter: Maximum number of iterations.  Default value: 300.
% * opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * opts.verbose: Flag for print out output, logical 0 or 1. Default value: 0.

%% Description
% The sg_min function has some optional input parameters that control the
% iteration process.  These parameters include maximum number of iteration,
% tolerance parameters for convergence of the objective function F and the 
% derivative of the objective function dF, and the print out of the iteration
% process.  The user can set one or all of parameters, if not, default
% values will be used.
% 

function opts = make_opts(opts)


if isfield(opts,'maxIter')
    if (opts.maxIter<1)
        opts.maxIter=300;
    end
else
    opts.maxIter=300;
end

if ~isfield(opts,'ftol')
    opts.ftol=1e-10;
end

if ~isfield(opts,'gradtol')
    opts.gradtol=1e-7;
end

if isfield(opts,'verbose')
    if (opts.verbose~=1)
        opts.verbose=0;
    end
else
    opts.verbose=0;
end

%% make_opts
% Make optional input parameters for running the sg_min package.

%% Syntax
%         Opts = make_opts(Opts)
%
%% Input 
%
% *Opts*: A list containing optional input parameter for sg_min.m specified
% by users.  One or several (even all) fields could be empty.
%
% * Opts.maxIter: Maximum number of iterations.
% * Opts.ftol: Tolerance parameter for F.  
% * Opts.gradtol: Tolerance parameter for dF.  
% * Opts.verbose: Flag for print out output, logical 0 or 1. 
%
%% Output:
% 
% *Opts*: A list containing optional input parameter for sg_min.m, specified by
% users or the default values are used.
%
% * Opts.maxIter: Maximum number of iterations.  Default value: 300.
% * Opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * Opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * Opts.verbose: Flag for print out output, logical 0 or 1. Default value: 0.

%% Description
% The sg_min function has some optional input parameters that control the
% iteration process.  These parameters include maximum number of iteration,
% tolerance parameters for convergence of the objective function F and the 
% derivative of the objective function dF, and the print out of the iteration
% process.  The user can set one or all of parameters, if not, default
% values will be used.
% 

function Opts = make_opts(Opts)


if isfield(Opts, 'maxIter')
    if (Opts.maxIter < 1)
        Opts.maxIter = 300;
    end
else
    Opts.maxIter = 300;
end

if ~isfield(Opts, 'ftol')
    Opts.ftol = 1e-10;
end

if ~isfield(Opts, 'gradtol')
    Opts.gradtol = 1e-7;
end

if isfield(Opts, 'verbose')
    if (Opts.verbose ~= 1)
        Opts.verbose = 0;
    end
else
    Opts.verbose = 0;
end

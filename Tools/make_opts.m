function opts = make_opts(opts)

% Options for Envelope method
%
% Notice:
% If one or several (even all) fields are empty, make_opts shall assign the
% default settings.
%
% If some fields of opts have been defined, make_opts shall check the fields
% for possible errors.
%
%
% Table of Options.  * * indicates default value.
%
%% FIELD            DESCRIPTION
%
% .init             .init specifies the initial r * u matrix.  
%                       *0* => no intial value defined.
%
%
% .maxIter          Maximum number of iterations.
%                       *300*
%
% .ftol             Tolerance parameter for F.
%                       *1e-10*
%
% .gradtol            Tolerance parameter for dF.
%                       *1e-7*
%
%
%
% .verbose          Flag for print out output.
%
%
%
%

% if ~isfield(opts,'init')
%     opts.init=[];
% end



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

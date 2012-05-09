function ddf = ddF(dF,Y,H)
% ddf = DDF(dF,Y,H)
%
% This function computes the second derivative of F, that is,
% ddf = d/dx dF(Y+x*H).
%
% dF is a handle to the function derivative.
% Y is expected to satisfy Y'*Y = I
% H is expected to be the same size as Y
% ddf will be the same size as Y
%
% See SG_MIN documentation for further details.
% =============================================================

ep = 1e-6;
n = norm(H,'fro');
ddf = (dF(Y+ep*H/n)-dF(Y-ep*H/n))/(2*ep/n);

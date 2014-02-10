% SG_MIN
% MATLAB Version 8.0 (R2012b) 12-Feb-2014
%
% Files
%   connection - Produces the christoffel symbol, C, for either the canonical or euclidean connections at the point Y
%   dimension  - Correctly counts the dimension, dim (stiefel dimension minus the grassmann degrees of freedom)
%   dtangent   - Computes the differential of the tangent map.  That is T = d/dt tangent(params,Y+t*dY,H)
%   ip         - Computes the inner produce of H1,H2 which are tangents at the stiefel point Y
%   move       - Moves the point Yi and direction Hi along a geodesic of length t in the metric to a point Yo and a direction Ho
%   nosym      - Orthogonal projection to remove the block diagonal components of the tangent vector H corresponding to the block diagonal right symmetries of F(Y)
%   tangent    - Produces tangent H from unconstrained differential D. H will be the unique tangent vector such that for any tangent W, real(trace(D'*W)) = ip(Y,W,H).  This is not always an orthogonal projection onto the tangent space

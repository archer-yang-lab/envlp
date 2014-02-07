% SG_MIN_REV
%
% Files
%   clamp           - reduces small roundoff errors in Y and D which can
%   ddF             - ddf = DDF(dF,Y,H)
%   dgrad           - Computes the tangent vector, W, which results from applying the
%   grad            - computes the gradient of the energy F at the point Y.
%   gradline        - the square magnitude of the gradient of the objective function,
%   invdgrad_CG     - Inverts the operator dgrad.  i.e. solves for H satisfying
%   invdgrad_MINRES - Inverts the operator dgrad.  i.e. solves for H satisfying
%   Kinetic         - 
%   make_dF         - make_dF
%   make_dFline     - DFLINE the derivative of the objective function, F, 
%   make_F          - make_F
%   make_Fline      - FLINE	the objective function, F, along the geodesic passing through Y
%   partition       - makes a guess at the partition of the function F
%   sg_dog          - SG_DOG(Y)	Optimize the objective function, F(Y) over all
%   sg_frcg         - SG_FRCG(Y)	Optimize the objective function, F(Y) over all
%   sg_min          - Stiefel/Grassmann minimization meant for template use.
%   sg_newton       - SG_NEWTON(Y)	Optimize the objective function, F(Y) over all
%   sg_prcg         - SG_FRCG(Y)	Optimize the objective function, F(Y) over all

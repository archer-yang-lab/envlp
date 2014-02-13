% SG_MIN_REV
% MATLAB Version 8.0 (R2012b) 12-Feb-2014
%
% Files
%   clamp           - Reduces small roundoff errors in Y and D which can accumulate during covariant operations that detroy the orthogonality conditions
%   ddF             - This function computes the second derivative of F, that is, ddf = d/dx dF(Y+x*H)
%   dgrad           - Computes the tangent vector, W, which results from applying the geometrically correct hessian at the stiefel point,Y, to the tangent vector H
%   grad            - Computes the gradient of the energy F at the point Y.
%   gradline        - GRADLINE the square magnitude of the gradient of the objective function, F, along the geodesic passing through Y in the direction H a distance t
%   invdgrad_CG     - Inverts the operator dgrad.  i.e. solves for H satisfying dl*H+dgrad(Y,H) = W.  The parameter, dl, is used for dogleg steps, and defaults to 0.  Uses a CG algorithm with a tolerance of gep*norm(W), or a tolerance of tol if given, if tol is given
%   invdgrad_MINRES - Inverts the operator dgrad.  i.e. solves for H satisfying dl*H+dgrad(Y,H) = W.  The parameter, dl, is used for dogleg steps, and defaults to 0.  Uses a MINRES algorithm with a tolerance of gep*norm(W), or a tolerance of tol if given, if tol is given
%   Kinetic         - Generates a square tridiagonal matrix, with the main diagonal elements all being 2, subdiagonal and superdiagonal elements all being -1
%   make_dF         - Generic function to generate the derivative function of the objective function F
%   make_dFline     - DFLINE the derivative of the objective function, F, along the geodesic passing through Y in the direction H a distance t
%   make_F          - Generic function to generate the objective function F
%   make_Fline      - FLINE	the objective function, F, along the geodesic passing through Y in the direction H a distance t
%   partition       - makes a guess at the partition of the function F based on the properties of dF at a randomly selected point Y
%   sg_dog          - Optimize the objective function, F(Y) over all Y such that Y'*Y = I.  Employs a local iterative search with initial point Y and terminates if the magnitude of the gradient falls to gradtol*(initial gradient magnitude) or if the relative decrease in F after some iteration is less than ftol
%   sg_frcg         - Optimize the objective function, F(Y) over all Y such that Y'*Y = I.  Fletcher-Reeves CG iterative search with initial point Y and terminates if the magnitude of the gradient falls to gradtol*(initial gradient magnitude) or if the relative decrease in F after some iteration is less than ftol
%   sg_min          - Stiefel/Grassmann minimization meant for template use.
%   sg_newton       - Optimize the objective function, F(Y) over all Y such that Y'*Y = I.  Employs a local iterative search with initial point Y and terminates if the magnitude of the gradient falls to gradtol*(initial gradient magnitude) or if the relative decrease in F after some iteration is less than ftol
%   sg_prcg         - Optimize the objective function, F(Y) over all Y such that Y'*Y = I.  Polak-Ribiere CG iterative search with initial point Y and terminates if the magnitude of the gradient falls to gradtol*(initial gradient magnitude) or if the relative decrease in F after some iteration is less than ftol
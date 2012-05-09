function Fline = make_Fline(F)
% FLINE	the objective function, F, along the geodesic passing through Y
%	in the direction H a distance t.
%
%	f= FLINE(t,Y,H)
%	Y is expected to satisfy Y'*Y = I.
%	H is expected to satisfy H = tangent(Y,H0)
%	
% role	high level algorithm, basic window dressing for the fmin function
Fline = @tmp;
    function fval = tmp(t,Y,H)
    	fval = F(move(Y,H,t));
    end
end

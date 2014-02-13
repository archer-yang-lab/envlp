function [fn,Yn]= sg_newton(F,dF,Y)
% SG_NEWTON(Y)	Optimize the objective function, F(Y) over all
%	Y such that Y'*Y = I.  Employs a local iterative search with
%	initial point Y and terminates if the magnitude of the gradient 
%	falls to gradtol*(initial gradient magnitude) or if the relative 
%	decrease in F after some iteration is less than ftol.
%
%	[fn,Yn]= SG_NEWTON(Y)
%	Y is expected to satisfy Y'*Y = I.
%	Yn will satisfy Yn'*Yn = I.
% role	high level algorithm, Newton's Methods
    global SGParameters;
    gradtol = SGParameters.gradtol;
    ftol = SGParameters.ftol;
    maxiter = SGParameters.maxiter;

    if (SGParameters.verbose)
	global SGdata;
	SGdata=[];
    end
    g = grad(Y); mag = sqrt(ip(Y,g,g));
    geps = mag*gradtol;
    f = F(Y);
    feps = ftol;

    N = 0; oldf = 2*f; oldmag = mag;
    if (SGParameters.verbose)
	SGdata = [];
	disp(sprintf('%s\t%s\t\t%s\t\t%s','iter','grad','F(Y)','step type'));
	SGdata(N+1,:) = [N mag f];
	disp(sprintf('%d\t%e\t%e\t%9d\t%s',N,mag,f,'none'));
    end

    while ((mag>geps) || (abs(oldf/f-1)>feps)) && (N < maxiter)
	N= N+1;

	gradsat = (mag<=geps);
	fsat = (abs(oldf/f-1)<=feps);
	if (fsat) fun = 'gradline'; else fun = 'Fline'; end

	sdr = -g; 
	gsdr = ip(Y,sdr,g);
	sa = -f/gsdr;
	sa = fminbnd(fun,-2*sa,2*sa,[],Y,sdr);
	sdrHsdr = ip(Y,sdr,dgrad(Y,sdr));
	sb = -gsdr/abs(sdrHsdr);
	sb = fminbnd(fun,-2*sb,2*sb,[],Y,sdr);
	Ysa = move(Y,sdr,sa);
	Ysb = move(Y,sdr,sb);

% note: the MINRES algorithm is more reliable of stiff problems than the CG
% algorithm.
	ndr = invdgrad_CG(Y,-g,gradtol*oldmag); 
	gndr = ip(Y,ndr,g);
	na = -f/gndr;
	na = fminbnd(fun,-2*na,2*na,[],Y,ndr);
	ndrHndr = ip(Y,ndr,dgrad(Y,ndr));
	nb = -gndr/abs(ndrHndr);
	nb = fminbnd(fun,-2*nb,2*nb,[],Y,ndr);
	Yna = move(Y,ndr,na);
	Ynb = move(Y,ndr,nb);

	if (fsat)
		gsa = grad(Ysa);
		gsb = grad(Ysb);
		gna = grad(Yna);
		gnb = grad(Ynb);
		magsa = sqrt(ip(Ysa,gsa,gsa)); magsb = sqrt(ip(Ysb,gsb,gsb));
		magna = sqrt(ip(Yna,gna,gna)); magnb = sqrt(ip(Ynb,gnb,gnb));
		if (min(magsa,magsb)<min(magna,magnb))
			if (SGParameters.verbose) steptype='steepest step'; end
			if (magsa<magsb) mag = magsa; Y = Ysa; g = gsa;
			else mag = magsb; Y = Ysb; g = gsb; end
		else
			if (SGParameters.verbose) steptype='Newton step'; end
			if (magna<magnb) mag = magna; Y = Yna; g = gna;
			else mag = magnb; Y = Ynb; g = gnb; end
		end
		newf = F(Y);
	else
		fsa = F(Ysa);
		fsb = F(Ysb);
		fna = F(Yna);
		fnb = F(Ynb);
		if (min(fsa,fsb)<min(fna,fnb))
			if (SGParameters.verbose) steptype='steepest step'; end
			if (fsa<fsb) newf = fsa; Y = Ysa; 
			else newf= fsb; Y = Ysb; end
		else
			if (SGParameters.verbose) steptype='Newton step'; end
			if (fna<fnb) newf = fna; Y = Yna;
			else newf = fnb; Y = Ynb; end
		end
		g = grad(Y); mag = sqrt(ip(Y,g,g));
	end
	
	oldf = f; f = newf;
	mag=sqrt(ip(Y,g,g));
	if (SGParameters.verbose)
		SGdata(N+1,:) = [N mag f];
		disp(sprintf('%d\t%e\t%e\t%9d\t%s',N,mag,f,steptype));
	end
    end

    fn = f;
    Yn = Y;
    if (N==maxiter)
        disp('WARNING: reached maximum number of iterations without convergence for specified tolerances');
        disp(strcat('Magnitude of gradient is :',num2str(mag)));
        disp(strcat('Increment of F is :',num2str(f-oldf)));
    end


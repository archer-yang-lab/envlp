function [fn,Yn]= sg_frcg(F,dF,Y)
% SG_FRCG(Y)	Optimize the objective function, F(Y) over all
%	Y such that Y'*Y = I.  Fletcher-Reeves CG iterative search with
%	initial point Y and terminates if the magnitude of the gradient 
%	falls to gradtol*(initial gradient magnitude) or if the relative 
%	decrease in F after some iteration is less than ftol.
%
%	[fn,Yn]= SG_FRCG(Y)
%	Y is expected to satisfy Y'*Y = I.
%	Yn will satisfy Yn'*Yn = I.
% role	high level algorithm, Fletcher-Reeves Method
    global SGParameters;
    gradtol = SGParameters.gradtol;
    ftol = SGParameters.ftol;
    dim = SGParameters.dimension;
    maxiter = SGParameters.maxiter;
    Fline = make_Fline(F);
    dFline = make_dFline(dF);
    
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
	disp(sprintf('%s\t%s\t\t%s\t\t%s','iter','grad','F(Y)'));
	SGdata(N+1,:) = [N mag f];
	disp(sprintf('%d\t%e\t%e\t%9d',N,mag,f));
    end

    reset=1;
    while ((mag>geps) || (abs(oldf/f-1)>feps) || reset) && (N < maxiter)
	N= N+1;

	gradsat = (mag<=geps);
	fsat = (abs(oldf/f-1)<=feps);

	dr = -g; 
	if (~reset)
		Hessolddr = dgrad(Y,olddr);
		alpha = ip(Y,dr,Hessolddr)/ip(Y,olddr,Hessolddr);
		dr = dr-alpha*olddr;
	else reset=0; end

	gdr = ip(Y,dr,-g); dr = dr*sign(gdr); gdr=abs(gdr);
	Hessdr = dgrad(Y,dr); drHdr = ip(Y,dr,Hessdr);
	cga = abs(gdr/drHdr);
	if (fsat)
		cgb = fzero(dFline,[0,2*cga],[],Y,dr);
	else
		cgb = fminbnd(Fline,-cga,cga,[],Y,dr);
	end

	[Ya,dra] = move(Y,dr,cga);
	[Yb,drb] = move(Y,dr,cgb);
	if (fsat)
		ga = grad(Ya);
		gb = grad(Yb);
		maga = sqrt(ip(Ya,ga,ga)); magb = sqrt(ip(Yb,gb,gb));
		if (maga<magb) mag = maga; Y = Ya; g = ga; olddr=dra;
		else mag = magb; Y = Yb; g = gb; olddr = drb; end
		newf = F(Y);
	else
		fa = F(Ya);
		fb = F(Yb);
		if (fa<fb) newf = fa; Y = Ya; olddr = dra;
		else newf= fb; Y = Yb; olddr = drb; end
		g = grad(Y); mag = sqrt(ip(Y,g,g));
	end
 
	oldf= f; f = newf;
	if (SGParameters.verbose)
		SGdata(N+1,:) = [N mag f];
		disp(sprintf('%d\t%e\t%e',N,mag,f));
	end
	if (rem(N,dim)==0) reset=1; end
    end
    fn = f;
    Yn = Y;

    if (N==maxiter)
        disp('WARNING: reached maximum number of iterations without convergence for specified tolerances');
        disp(strcat('Magnitude of gradient is :',num2str(mag)));
        disp(strcat('Increment of F is :',num2str(f-oldf)));
    end


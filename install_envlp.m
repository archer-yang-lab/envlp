function  install_envlp

%INSTALL_ENVLP -- Add various envelope utilities to the MATLAB path
%
%  * SG_MIN Ver 2.4.3
%  * SG_MIN Ver 2.4.3 (Modified Components)
%  * env: Envelope
%  * envmean: Envelope for Estimating the Multivariate Mean
%  * envseq: Envelope Using Sequential Algorithm
%  * xenv: Envelope in the Predictor Space
%  * xenvpls: Envelope in the Predictor Space Using Partial Least Squares Algorithm
%  * henv: Heteroscedastic Envelope
%  * ienv: Inner Envelope
%  * penv: Partial Envelope
%  * senv: Scaled Envelope
%  * Post-pocessing Tools
%  * Auxiliary Functions
%  * Example Data
%  * Documentation
%
% See also ADDPATH, PERCLEARN_INSTALL, PATH2RC, REHASH.

% Original coding by Zhihua Su and Yi Yang, University of Minnesota
% $Revision: 1.0.0 $  $Date: 2012-05-03 $
%
% Part of the envlp toolbox version 1.0.0 for MATLAB version 5 and up.
% <http://code.google.com/p/envlp/ http://code.google.com/p/envlp/>
% Copyright (c) Dennis Cook, Zhihua Su, Yi Yang, 2012
% Please read the LICENSE and NO WARRANTY statement in ./envlp_license.m

%-- Print herald
more off ;
home ;


fprintf('envlp -- A MATLAB Toolbox for Envelope Models ver. 1.0.0.\n\n') ;
fprintf('Released DD MMM YYYY.\n\n') ;

fprintf('Copyright (c) Dennis Cook, Zhihua Su, Yi Yang, 2012.\n\n') ;
fprintf('<a href="http://code.google.com/p/envlp/">http://code.google.com/p/envlp/</a>.\n\n') ;
fprintf('This software is freely available and freely redistributable,\n') ;
fprintf('according to the conditions of the GNU General Public License.\n') ;
fprintf('The full text of the GNU General Public License is available at\n') ;
fprintf('the Free Software Foundation website (<a href="http://www.fsf.org">http://www.fsf.org</a>) and is\n') ;
fprintf('reproduced in envlp_license.m. In brief, the GNU GPL provisions are:\n\n') ;

fprintf('You may not distribute the software, in whole or in part, in\n') ;
fprintf('conjunction with proprietary code. That means you ONLY have my\n') ;
fprintf('permission to distribute a program that uses my code IF you also\n') ;
fprintf('make freely available (under the terms of the GNU GPL) the source\n') ;
fprintf('code for your whole project. You may not pass on the software to\n') ;
fprintf('another party in its current form or any altered, embellished or\n') ;
fprintf('reduced form, without acknowledging the author and including a\n') ;
fprintf('copy of the GNU General Public License.\n\n') ;

fprintf('The software is distributed in the hope that it will be useful,\n') ;
fprintf('but WITHOUT ANY WARRANTY; without even the implied warranty of\n') ;
fprintf('MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n') ;
fprintf('See the GNU General Public License for more details.\n\n') ;

reply = lower(input('Do you agree with this license agreement? y/n (y): ','s')) ;
if (isempty(reply))  reply = 'y' ; end
if (reply(1)~='y')
    fprintf('\nOK. No changes were made to the path. Good bye!\n\n') ;
    return
end
fprintf('\n#############################################################################\n') ;
fprintf('\n## NOTE: If you have previously installed an older version of this toolbox ##\n') ;
fprintf('\n## please have it removed from MATLAB search path first before installing  ##\n') ;
fprintf('\n## this version.                                                           ##\n') ;
fprintf('\n#############################################################################\n') ;
fprintf('\nThis script will install the following utilities to the Matlab path.\n') ;
fprintf(' * SG_MIN Ver 2.4.3\n') ;
fprintf(' * SG_MIN Ver 2.4.3 (Modified Components)\n') ;
fprintf(' * env: Envelope\n') ;
fprintf(' * envmean: Envelope for Estimating the Multivariate Mean\n');
fprintf(' * envseq: Envelope Using Sequential Algorithm\n');
fprintf(' * xenv: Envelope in the Predictor Space\n') ;
fprintf(' * xenvpls: Envelope in the Predictor Space Using Partial Least Squares Algorithm\n');
fprintf(' * henv: Heteroscedastic Envelope\n') ;
fprintf(' * ienv: Inner Envelope\n') ;
fprintf(' * penv: Partial Envelope\n') ;
fprintf(' * senv: Scaled Envelope\n') ;
fprintf(' * Auxiliary Tools\n') ;
fprintf(' * Example Data\n') ;
fprintf(' * Documentation\n') ;

fprintf('\nSee "/doc/tutorial.m" for a quick start; See "README.txt" for more help."\n') ;


reply = lower(input('Do you want to proceed? y/n (y): ','s')) ;
if (isempty(reply))  reply = 'y' ; end
if (reply(1)~='y')
    fprintf('\nOK. No changes were made to the path. Good bye!\n\n') ;
    return
end

%-- Path name -- edit this if necessary
envlp_path = fileparts(which('install_envlp')) ;  % e.g. '/work/envlp'

s = fullfile(envlp_path,'/src/sg_min') ;
addpath(s) ;
fprintf('\nSG_MIN Ver 2.4.3 added to the path: %s\n\n',s) ;

s = fullfile(envlp_path,'/src/sg_min_rev') ;
addpath(s) ;
fprintf('SG_MIN Ver 2.4.3 (Modified Components) added to the path: %s\n\n',s) ;

s = fullfile(envlp_path,'/src/env') ;
addpath(s) ;
fprintf('env added to the path: %s\n\n',s) ;

s = fullfile(envlp_path,'/src/envmean') ;
addpath(s) ;
fprintf('envmean added to the path: %s\n\n',s) ;

s = fullfile(envlp_path,'/src/envseq') ;
addpath(s) ;
fprintf('envseq added to the path: %s\n\n',s) ;

s = fullfile(envlp_path,'/src/xenv') ;
addpath(s) ;
fprintf('xenv added to the path: %s\n\n',s) ;

s = fullfile(envlp_path,'/src/xenvpls') ;
addpath(s) ;
fprintf('xenvpls added to the path: %s\n\n',s) ;


s = fullfile(envlp_path,'/src/henv') ;
addpath(s) ;
fprintf('henv added to the path: %s\n\n',s) ;

s = fullfile(envlp_path,'/src/ienv') ;
addpath(s) ;
fprintf('ienv added to the path: %s\n\n',s) ;

s = fullfile(envlp_path,'/src/penv') ;
addpath(s) ;
fprintf('penv added to the path: %s\n\n',s) ;

s = fullfile(envlp_path,'/src/senv') ;
addpath(s) ;
fprintf('senv added to the path: %s\n\n',s) ;

s = fullfile(envlp_path,'/src/tools') ;
addpath(s) ;
fprintf('Post-pocessing Tools added to the path: %s\n\n',s) ;

s = fullfile(envlp_path,'/src/auxiliary') ;
addpath(s) ;
fprintf('Auxiliary functions added to the path: %s\n\n',s) ;

s = fullfile(envlp_path,'data') ;
addpath(s) ;
fprintf('Example Data added to the path: %s\n\n',s) ;

s = fullfile(envlp_path,'doc') ;
addpath(s) ;
fprintf('Documentation added to the path: %s\n\n',s) ;

s = fullfile(envlp_path,'examples') ;
addpath(s) ;
fprintf('Examples files added to the path: %s\n\n',s) ;

%-- Save and finish
if (which('savepath'))
    notsaved = savepath ;  % Save the changes for future Matlab sessions.
else
    notsaved = path2rc ;   % deprecated version of SAVEPATH
end
if (notsaved)
    fprintf('\nThere were problems saving the path. ') ;
    fprintf('Manual intervention seems necessary.\n\n') ;
else
    rehash toolbox ;   % Make sure the new settings are in effect.
    rehash toolboxcache ;
    fprintf('The updated path saved to MATLAB''s pathdef.m file.\n\n') ;
end

fprintf('Updates may be available at <a href="http://code.google.com/p/envlp/">http://code.google.com/p/envlp/</a>\n') ;
fprintf('Thank you for your interest in this software.\n') ;

%---  Return no values
%%%%% End of file INSTALL_ENVELOPE.M

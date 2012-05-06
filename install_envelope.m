function  install_envelope

%INSTALL_ENVELOPE -- Add various envelope utilities to the MATLAB path
%
%  * SG_MIN Ver 2.4.1
%  * SG_MIN Ver 2.4.1 (Modified Components)
%  * Envelope Module
%  * Envelope Module for Predictor Reduction
%  * Heteroscedastic Envelope Module
%  * Inner Envelope Module
%  * Partial Envelope Module
%  * Scaled Envelope Module
%  * Auxiliary Tools
%  * Example Data
%  * Documentation
%
% See also ADDPATH, PERCLEARN_INSTALL, PATH2RC, REHASH.

% Original coding by Zhihua Su and Yi Yang, University of Minnesota
% $Revision: 1.0.0 $  $Date: 2012-05-03 $
%
% Part of the envelope toolbox version 1.1 for MATLAB version 5 and up.
% http://code.google.com/p/envlp/
% Copyright (c) Dennis Cook, Zhihua Su, Yi Yang, 2012
% Please read the LICENSE and NO WARRANTY statement in ./envelope_license.m

%-- Print herald
more off ;
home ;


fprintf('Matlab ENVELOPE Toolbox ver. 1.0.0, released DD MMM YYYY.\n\n') ;

fprintf('Copyright (c) Dennis Cook, Zhihua Su, Yi Yang, 2012.\n\n') ;
fprintf('http://code.google.com/p/envlp/.\n\n') ;
fprintf('This software is freely available and freely redistributable,\n') ;
fprintf('according to the conditions of the GNU General Public License.\n') ;
fprintf('The full text of the GNU General Public License is available at\n') ;
fprintf('the Free Software Foundation website (http://www.fsf.org) and is\n') ;
fprintf('reproduced in envelope_license.m. In brief, the GNU GPL provisions are:\n\n') ;

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

fprintf('\nThis script will install the following utilities to the Matlab path.\n') ;
fprintf(' * SG_MIN Ver 2.4.1\n') ;
fprintf(' * SG_MIN Ver 2.4.1 (Modified Components)\n') ;
fprintf(' * Envelope Module\n') ;
fprintf(' * Envelope Module for Predictor Reduction\n') ;
fprintf(' * Heteroscedastic Envelope Module\n') ;
fprintf(' * Inner Envelope Module\n') ;
fprintf(' * Partial Envelope Module\n') ;
fprintf(' * Scaled Envelope Module\n') ;
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
envelope_path = fileparts(which('install_envelope')) ;  % e.g. '/work/envelope'

s = fullfile(envelope_path,'/src/sg_min') ;
addpath(s) ;
fprintf('\nSG_MIN Ver 2.4.1 added to the path: %s\n\n',s) ;

s = fullfile(envelope_path,'/src/sg_min_rev') ;
addpath(s) ;
fprintf('SG_MIN Ver 2.4.1 (Modified Components) added to the path: %s\n\n',s) ;

s = fullfile(envelope_path,'/src/env') ;
addpath(s) ;
fprintf('Envelope Module added to the path: %s\n\n',s) ;

s = fullfile(envelope_path,'/src/xenv') ;
addpath(s) ;
fprintf('Envelope Module for Predictor Reduction added to the path: %s\n\n',s) ;

s = fullfile(envelope_path,'/src/henv') ;
addpath(s) ;
fprintf('Heteroscedastic Envelope Module added to the path: %s\n\n',s) ;

s = fullfile(envelope_path,'/src/ienv') ;
addpath(s) ;
fprintf('Inner Envelope Module added to the path: %s\n\n',s) ;

s = fullfile(envelope_path,'/src/penv') ;
addpath(s) ;
fprintf('Partial Envelope Module added to the path: %s\n\n',s) ;

s = fullfile(envelope_path,'/src/senv') ;
addpath(s) ;
fprintf('Scaled Envelope Module added to the path: %s\n\n',s) ;

s = fullfile(envelope_path,'/src/tools') ;
addpath(s) ;
fprintf('Auxiliary Tools added to the path: %s\n\n',s) ;

s = fullfile(envelope_path,'data') ;
addpath(s) ;
fprintf('Example Data added to the path: %s\n\n',s) ;

s = fullfile(envelope_path,'doc') ;
addpath(s) ;
fprintf('Documentation added to the path: %s\n\n',s) ;

s = fullfile(envelope_path,'examples') ;
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

fprintf('Updates may be available at http://code.google.com/p/envlp/ \n') ;
fprintf('Thank you for your interest in this software.\n') ;

%---  Return no values
%%%%% End of file INSTALL_ENVELOPE.M

% This example demonstrates the usage of xenvpls module.

load wheatprotein.txt  % Load data
X = wheatprotein(:, 1 : 6); % Assign predictors
Y = wheatprotein(:, 7); % Assign response


% Dimension selection
m = 5;
u = mfoldcv(X, Y, m, 'xenvpls') % Select u using 5-fold cross validation

% % u =
% % 
% %      0


% Fit the model
ModelOutput = xenvpls(X, Y, u)

% % ModelOutput = 
% % 
% %       beta: [6x1 double]
% %       SigX: [6x6 double]
% %      Gamma: []
% %     Gamma0: [6x6 double]
% %        eta: []
% %      Omega: []
% %     Omega0: [6x6 double]
% %         mu: 10.0538
% %     sigYcX: 2.1084
% %         np: 23
% %          n: 50


ModelOutput.beta  % Print estimated regression coefficients

% % ans =
% % 
% %      0
% %      0
% %      0
% %      0
% %      0
% %      0
   
   
ModelOutput.sigYcX  % Check error variance

% % ans =
% % 
% %     2.1084


% Inference tools
B = 100;
bootse = bootstrapse(X, Y, u, B, 'xenvpls') % Compute bootstrap standard errors for elements in \beta

% % bootse =
% % 
% %      0
% %      0
% %      0
% %      0
% %      0
% %      0



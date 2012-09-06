% This example demonstrates the usage of xenvpls module.

load wheatprotein.txt  % Load data
X = wheatprotein(:, 1 : 6); % Assign predictors
Y = wheatprotein(:, 7); % Assign response


% Dimension selection
m = 5;
u = mfoldcv(X, Y, m, 'xenvpls') % Select u using 5-fold cross validation

% % u =
% % 
% %      6


% Fit the model
ModelOutput = xenvpls(X, Y, u)

% % ModelOutput = 
% % 
% %       beta: [6x1 double]
% %       SigX: [6x6 double]
% %      Gamma: [6x6 double]
% %     Gamma0: []
% %        eta: [6x1 double]
% %      Omega: [6x6 double]
% %     Omega0: []
% %         mu: 24.5781
% %     sigYcX: 0.0321
% %         np: 29
% %          n: 50


ModelOutput.beta  % Print estimated regression coefficients

% % ans =
% % 
% %    -0.0416
% %    -0.0490
% %     0.3368
% %    -0.1981
% %     0.0020
% %    -0.0480
   
   
ModelOutput.sigYcX  % Check error variance

% % ans =
% % 
% %     0.0321


% Inference tools
B = 100;
bootse = bootstrapse(X, Y, u, B, 'xenvpls') % Compute bootstrap standard errors for elements in \beta

% % bootse =
% % 
% %     0.0255
% %     0.0336
% %     0.0349
% %     0.0179
% %     0.0021
% %     0.0095



% This example demonstrates the usage of envseq module.

load Rohwer  % Load data
X = Rohwer(:, 4 : 5); % Assign predictors
Y = Rohwer(:, 1 : 3); % Assign responses


% Dimension selection
m = 5;
u = mfoldcv(X, Y, m, 'envseq') % Select u using 5-fold cross validation

% % u =
% % 
% %      1


% Fit the model
ModelOutput = envseq(X, Y, u)

% % ModelOutput = 
% % 
% %         beta: [3x2 double]
% %        Sigma: [3x3 double]
% %        Gamma: [3x1 double]
% %       Gamma0: [3x2 double]
% %          eta: [2.5708 1.1966]
% %        Omega: 752.8146
% %       Omega0: [2x2 double]
% %        alpha: [3x1 double]
% %     paramNum: 11
% %            n: 69


ModelOutput.beta  % Print estimated regression coefficients

% % ans =
% % 
% %     2.0952    0.9752
% %     1.4774    0.6877
% %     0.1911    0.0889
   
   
ModelOutput.Sigma  % Check error covariance matrix

% % ans =
% % 
% %   587.2575  229.9647   37.2716
% %   229.9647  421.1372   42.9256
% %    37.2716   42.9256   12.2696


% Inference tools
B = 100;
bootse = bootstrapse(X, Y, u, B, 'envseq') % Compute bootstrap standard errors for elements in \beta

% % bootse =
% % 
% %     0.8834    0.6340
% %     0.5204    0.4381
% %     0.0885    0.0616



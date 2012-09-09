% This example demonstrates the usage of xenvpls module.

load VocabGrowth % Load data
X = VocabGrowth(:, 1 : 3); % Assign predictors
Y = VocabGrowth(:, 4); % Assign response


% Dimension selection
m = 5;
u = mfoldcv(X, Y, m, 'xenvpls') % Select u using 5-fold cross validation

% % u =
% % 
% %      1


% Fit the model
ModelOutput = xenvpls(X, Y, u)

% % ModelOutput = 
% % 
% %         beta: [3x1 double]
% %         SigX: [3x3 double]
% %        Gamma: [3x1 double]
% %       Gamma0: [3x2 double]
% %          eta: 5.2899
% %        Omega: 10.9286
% %       Omega0: [2x2 double]
% %           mu: 1.5683
% %       sigYcX: 1.0934
% %     paramNum: 9
% %            n: 64


ModelOutput.beta  % Print estimated regression coefficients

% % ans =
% % 
% %     0.2573
% %     0.2741
% %     0.3049
   
   
ModelOutput.sigYcX  % Check error variance

% % ans =
% % 
% %     1.0934


% Inference tools
B = 100;
bootse = bootstrapse(X, Y, u, B, 'xenvpls') % Compute bootstrap standard errors for elements in \beta

% % bootse =
% % 
% %     0.0219
% %     0.0242
% %     0.0259



% This example demonstrates the usage of envseq module.

load wheatprotein.txt  % Load data
X = wheatprotein(:, 8); % Assign predictor
Y = wheatprotein(:, 1 : 6); % Assign responses


% Dimension selection
m = 5;
u = mfoldcv(X, Y, m, 'envseq') % Select u using 5-fold cross validation

% % u =
% % 
% %      0


% Fit the model
ModelOutput = envseq(X, Y, u)

% % ModelOutput = 
% % 
% %       beta: [6x1 double]
% %      Sigma: [6x6 double]
% %      Gamma: []
% %     Gamma0: [6x6 double]
% %        eta: []
% %      Omega: []
% %     Omega0: [6x6 double]
% %      alpha: [6x1 double]
% %   paramNum: 27
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
   
   
ModelOutput.Sigma  % Check error covariance matrix

% % ans =
% % 
% %    1.0e+03 *
% % 
% %   Columns 1 through 4
% % 
% %     1.1932    0.9825    1.0568    1.1508
% %     0.9825    0.8222    0.8808    0.9391
% %     1.0568    0.8808    0.9451    1.0125
% %     1.1508    0.9391    1.0125    1.1239
% %     1.5411    1.2700    1.3725    1.5217
% %     0.6340    0.5279    0.5652    0.6180
% % 
% %   Columns 5 through 6
% % 
% %     1.5411    0.6340
% %     1.2700    0.5279
% %     1.3725    0.5652
% %     1.5217    0.6180
% %     2.3255    0.8365
% %     0.8365    0.3619


% Inference tools
B = 100;
bootse = bootstrapse(X, Y, u, B, 'envseq') % Compute bootstrap standard errors for elements in \beta

% % bootse =
% % 
% %      0
% %      0
% %      0
% %      0
% %      0
% %      0



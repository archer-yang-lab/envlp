% This example demonstrates the usage of xenv module.
% The data is used as an example in Cook et al. (2012).

load wheatprotein.txt  % Load data
X = wheatprotein(:, 1 : 6); % Assign predictors
Y = wheatprotein(:, 7); % Assign response


% Dimension selection
u1 = modelselectaic(X, Y, 'xenv') % Select u using AIC

% % u1 =
% % 
% %      4


u2 = modelselectbic(X, Y, 'xenv') % Select u using BIC

% % u2 =
% % 
% %      4

     
alpha = 0.01;
u3 = modelselectlrt(X, Y, alpha, 'xenv') % Select u using LRT

% % u3 =
% % 
% %      4
     

% Fit the model
ModelOutput = xenv(X, Y, u3)

% % ModelOutput = 
% % 
% %          beta: [6x1 double]
% %          SigX: [6x6 double]
% %         Gamma: [6x4 double]
% %        Gamma0: [6x2 double]
% %           eta: [4x1 double]
% %         Omega: [4x4 double]
% %        Omega0: [2x2 double]
% %            mu: 24.8863
% %        sigYcX: 0.0321
% %             l: -865.6407
% %     covMatrix: [6x6 double]
% %       asyXenv: [6x1 double]
% %         ratio: [6x1 double]
% %            np: 27
% %             n: 50


ModelOutput.beta % Print estimated regression coefficients

% % ans =
% % 
% %    -0.0443
% %    -0.0481
% %     0.3377
% %    -0.1963
% %     0.0019
% %    -0.0487
   
   
ModelOutput.sigYcX % Check the error variance

% % ans =
% % 
% %    0.0321
    

eig(ModelOutput.Omega) % Check the eigenvalues of \Omega

% % ans =
% % 
% %   208.4944
% %     0.2664
% %    25.9938
% %    20.0812
    
    
eig(ModelOutput.Omega0) % Check the eigenvalues of \Omega_0

% % ans =
% % 
% %    1.0e+03 *
% % 
% %     0.0004
% %     6.5166

    
% Inference tools
B = 100;
bootse = bootstrapse(X, Y, u3, B, 'xenv') % Compute bootstrap standard errors for elements in \beta

% % bootse =
% % 
% %     0.0265
% %     0.0376
% %     0.0422
% %     0.0177
% %     0.0025
% %     0.0088


Xnew = mean(X)';
PredictOutput = prediction(ModelOutput, Xnew, 'prediction', 'xenv') % Prediction based on a new observation

% % PredictOutput = 
% % 
% %         value: 10.0538
% %     covMatrix: 16.8778
% %            SE: 4.1083


PredictOutput.value % Check predicted value

% % ans =
% % 
% %   10.0538


PredictOutput.SE % Check the standard errors of each elements in the predicted value

% % ans =
% % 
% %    4.1083


TestOutput = testcoefficient(ModelOutput, 'xenv'); % Test if \beta = 0

% %  Test Hypothesis     Chisq Statistic    DF     P-value
% % -------------------------------------------------------
% % L * beta * R = A        3233.053        6      0.0000
% % -------------------------------------------------------


TestOutput

% % TestOutput = 
% % 
% %     chisqStatistic: 3.2331e+03
% %                 df: 6
% %             pValue: 0
% %          covMatrix: [6x6 double]


TestInput.L = ones(1, 6); % Define L
TestInput.A = 0; % Define A
TestOutput = testcoefficient(ModelOutput, 'xenv', TestInput); % Test if L * \beta = A

% %  Test Hypothesis     Chisq Statistic    DF     P-value
% % -------------------------------------------------------
% % L * beta * R = A           0.681        1      0.4091
% % -------------------------------------------------------



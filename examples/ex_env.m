% This example demonstrates the usage of env module.
% The data is used as an example in Cook et al. (2010).

load wheatprotein.txt  % Load data
X = wheatprotein(:, 8); % Assign predictor
Y = wheatprotein(:, 1 : 6); % Assign responses


% Dimension selection
u1 = modelselectaic(X, Y, 'env') % Select u using AIC

% % u1 =
% % 
% %      1


u2 = modelselectbic(X, Y, 'env') % Select u using BIC

% % u2 =
% % 
% %      1

     
alpha = 0.01;
u3 = modelselectlrt(X, Y, alpha, 'env') % Select u using LRT

% % u3 =
% % 
% %      1
     

% Fit the model
ModelOutput = env(X, Y, u3)

% % ModelOutput = 
% % 
% %          beta: [6x1 double]
% %         Sigma: [6x6 double]
% %         Gamma: [6x1 double]
% %        Gamma0: [6x5 double]
% %           eta: 8.5647
% %         Omega: 7.8762
% %        Omega0: [5x5 double]
% %         alpha: [6x1 double]
% %             l: -850.7592
% %     covMatrix: [6x6 double]
% %        asyEnv: [6x1 double]
% %         ratio: [6x1 double]
% %      paramNum: 28
% %             n: 50


ModelOutput.beta % Print estimated regression coefficients

% % ans =
% % 
% %    -1.0644
% %     4.4730
% %     3.6839
% %    -5.9770
% %     0.6013
% %    -1.5986
   
   
ModelOutput.ratio % Check the relative efficiency of the envelope model versus the standard model

% % ans =
% % 
% %    28.0945
% %    18.4326
% %    23.6384
% %    16.3211
% %    65.8245
% %     6.4668
    

eig(ModelOutput.Omega) % Check the eigenvalues of \Omega

% % ans =
% % 
% %     7.8762
    
    
eig(ModelOutput.Omega0) % Check the eigenvalues of \Omega_0

% % ans =
% % 
% %    1.0e+03 *
% % 
% %     6.5166
% %     0.2083
% %     0.0201
% %     0.0004
% %     0.0003

    
% Inference tools
B = 100;
bootse = bootstrapse(X, Y, u3, B, 'env') % Compute bootstrap standard errors for elements in \beta

% % bootse =
% % 
% %     0.2955
% %     0.4082
% %     0.3568
% %     0.5454
% %     0.2228
% %     0.6600


Xnew = 1;
PredictOutput = prediction(ModelOutput, Xnew, 'prediction', 'env') % Prediction based on a new observation

% % PredictOutput = 
% % 
% %         value: [6x1 double]
% %     covMatrix: [6x6 double]
% %            SE: [6x1 double]


PredictOutput.value % Check predicted value

% % ans =
% % 
% %   473.6491
% %   131.9470
% %   254.8883
% %   374.8510
% %   381.5486
% %    -7.9273


PredictOutput.SE % Check the standard errors of each elements in the predicted value

% % ans =
% % 
% %    34.9178
% %    28.7313
% %    30.8796
% %    33.9056
% %    48.6949
% %    19.2628


TestOutput = testcoefficient(ModelOutput, 'env'); % Test if \beta = 0

% %  Test Hypothesis     Chisq Statistic    DF     P-value
% % -------------------------------------------------------
% % L * beta * R = A         116.230        6      0.0000
% % -------------------------------------------------------


TestOutput

% % TestOutput = 
% % 
% %     chisqStatistic: 116.2299
% %                 df: 6
% %             pValue: 0
% %          covMatrix: [6x6 double]


TestInput.L = ones(1, 6); % Define L
TestInput.A = 0; % Define A
TestOutput = testcoefficient(ModelOutput, 'env', TestInput); % Test if L * \beta = A

% %  Test Hypothesis     Chisq Statistic    DF     P-value
% % -------------------------------------------------------
% % L * beta * R = A           0.084        1      0.7724
% % -------------------------------------------------------


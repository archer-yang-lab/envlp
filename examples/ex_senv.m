% This example demonstrates the usage of senv module.
% The data is used as an example in Cook and Su (2012).

load('sales.txt')  % Load data
Y = sales(:, 4 : 7);  % Assign predictor
X = sales(:, 1 : 3);  % Assign responses


% Dimension selection
% With senv, we cannot apply LRT to select the dimension
u1 = modelselectaic(X, Y, 'senv') % Select u using AIC

% % u1 =
% % 
% %      4


u2 = modelselectbic(X, Y, 'senv') % Select u using BIC

% % u2 =
% % 
% %      2

     
% Fit the model
ModelOutput = senv(X, Y, u2)

% % ModelOutput = 
% % 
% %          beta: [4x3 double]
% %         Sigma: [4x4 double]
% %        Lambda: [4x4 double]
% %         Gamma: [4x2 double]
% %        Gamma0: [4x2 double]
% %           eta: [2x3 double]
% %         Omega: [2x2 double]
% %        Omega0: [2x2 double]
% %         alpha: [4x1 double]
% %            np: 23
% %             l: -386.1900
% %     covMatrix: [12x12 double]
% %       asySenv: [4x3 double]
% %         ratio: [4x3 double]
% %             n: 50


ModelOutput.beta % Print estimated regression coefficients

% % ans =
% % 
% %     0.2373   -0.0559    0.3030
% %    -0.0062    0.2656   -0.0229
% %     0.2651   -0.2197    0.3473
% %     0.3533    0.5755    0.4139
   
   
ModelOutput.ratio % Check the relative efficiency of the scaled envelope model versus the standard model

% % ans =
% % 
% %     3.3560    2.0912    2.8712
% %     2.5888    1.4747    2.2799
% %     1.7983    1.2281    1.6458
% %     2.1622    1.8652    2.0954
    

eig(ModelOutput.Omega) % Check the eigenvalues of \Omega

% % ans =
% % 
% %     0.4382
% %     1.1006
    
    
eig(ModelOutput.Omega0) % Check the eigenvalues of \Omega_0

% % ans =
% % 
% %     5.3040
% %    13.1733

    
% Inference tools
B = 100;
bootse = bootstrapse(X, Y, u2, B, 'senv') % Compute bootstrap standard errors for elements in \beta

% % bootse =
% % 
% %     0.1046    0.1209    0.1398
% %     0.0815    0.1042    0.1186
% %     0.0663    0.0846    0.0943
% %     0.1316    0.1252    0.1392


Xnew = mean(X)';
PredictOutput = prediction(ModelOutput, Xnew, 'estimation', 'senv') % Estimation of Y given the X of the observation

% % PredictOutput = 
% % 
% %         value: [4x1 double]
% %     covMatrix: [4x4 double]
% %            SE: [4x1 double]


PredictOutput.value % Check predicted value

% % ans =
% % 
% %    11.2200
% %    14.1800
% %    10.5600
% %    29.7600


PredictOutput.SE % Check the standard errors of each elements in the predicted value

% % ans =
% % 
% %     8.0897
% %     6.1029
% %     4.3341
% %     8.3830


TestInput.L = ones(1, 4); % Define L
TestInput.R = ones(3, 1); % Define R
TestInput.A = 1; % Define A
TestOutput = testcoefficient(ModelOutput, 'senv', TestInput); % Test if L * \beta * R = A

% %  Test Hypothesis     Chisq Statistic    DF     P-value
% % -------------------------------------------------------
% % L * beta * R = A         441.252        1      0.0000
% % -------------------------------------------------------


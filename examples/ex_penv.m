% This example demonstrates the usage of penv module.
% The data is used as an example in Su and Cook (2011).

load fiberpaper.dat  % Load data
Y = fiberpaper(:, 1 : 4);  % Assign responses
X.X1 = fiberpaper(:, 7);  % Assign the main predictor
X.X2 = fiberpaper(:, 5 : 6);  % Assign the covariates


% Dimension selection
u1 = modelselectaic(X, Y, 'penv') % Select u using AIC

% % u1 =
% % 
% %      3


u2 = modelselectbic(X, Y, 'penv') % Select u using BIC

% % u2 =
% % 
% %      1

     
alpha = 0.01;
u3 = modelselectlrt(X, Y, alpha, 'penv') % Select u using LRT

% % u3 =
% % 
% %      1
     

% Fit the model
ModelOutput = penv(X, Y, u3)

% % ModelOutput = 
% % 
% %         beta1: [4x1 double]
% %         beta2: [4x2 double]
% %         alpha: [4x1 double]
% %         Gamma: [4x1 double]
% %           eta: 0.0047
% %         Omega: 0.0149
% %        Omega0: [3x3 double]
% %         Sigma: [4x4 double]
% %             l: -35.6323
% %      paramNum: 23
% %     covMatrix: [12x12 double]
% %         asySE: [4x1 double]
% %         ratio: [4x1 double]
% %             n: 62


ModelOutput.beta1 % Print estimated coefficients for the main predictor(s)

% % ans =
% % 
% %    -0.0010
% %    -0.0030
% %     0.0034
% %     0.0008
   
   
ModelOutput.ratio % Check the relative efficiency of the partial envelope model versus the standard model

% % ans =
% % 
% %    65.9692
% %     6.8217
% %    10.4152
% %     9.6228
    

eig(ModelOutput.Omega) % Check the eigenvalues of \Omega

% % ans =
% % 
% %     0.0149
    
    
eig(ModelOutput.Omega0) % Check the eigenvalues of \Omega_0

% % ans =
% % 
% %     0.0050
% %     0.0999
% %     4.9819

    
% Inference tools
B = 100;
bootse = bootstrapse(X, Y, u3, B, 'penv') % Compute bootstrap standard errors for elements in \beta

% % bootse =
% % 
% %     0.0008
% %     0.0012
% %     0.0015
% %     0.0008


Xnew.X1 = X.X1(1, :)';
Xnew.X2 = X.X2(1, :)';
PredictOutput = prediction(ModelOutput, Xnew, 'estimation', 'penv') % Estimation of Y when X1 and X2 take the same value as the first sample 

% % PredictOutput = 
% % 
% %         value: [4x1 double]
% %     covMatrix: [4x4 double]
% %            SE: [4x1 double]


PredictOutput.value % Check predicted value

% % ans =
% % 
% %    21.1169
% %     7.1173
% %     5.3637
% %     0.8737


PredictOutput.SE % Check the standard errors of each elements in the predicted value

% % ans =
% % 
% %     1.4680
% %     0.4234
% %     0.7145
% %     0.3161


TestOutput = testcoefficient(ModelOutput, 'penv'); % Test if \beta1 = 0

% %  Test Hypothesis     Chisq Statistic    DF     P-value
% % -------------------------------------------------------
% % L * beta * R = A          12.604        4      0.0134
% % -------------------------------------------------------


TestOutput

% % TestOutput = 
% % 
% %     chisqStatistic: 12.6039
% %                 df: 4
% %             pValue: 0.0134
% %          covMatrix: [4x4 double]


TestInput.L = [1 2 3 4]; % Define L
TestInput.A = 0; % Define A
TestOutput = testcoefficient(ModelOutput, 'penv', TestInput); % Test if L * \beta = A

% %  Test Hypothesis     Chisq Statistic    DF     P-value
% % -------------------------------------------------------
% % L * beta * R = A           6.706        1      0.0096
% % -------------------------------------------------------


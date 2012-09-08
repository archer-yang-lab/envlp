% This example demonstrates the usage of ienv module.
% The data is used as an example in Su and Cook (2010).

load irisf.mat  % Load data


% Dimension selection
u1 = modelselectaic(X, Y, 'ienv') % Select u using AIC

% % u1 =
% % 
% %      1


u2 = modelselectbic(X, Y, 'ienv') % Select u using BIC

% % u2 =
% % 
% %      1

     
alpha = 0.01;
u3 = modelselectlrt(X, Y, alpha, 'ienv') % Select u using LRT

% % u3 =
% % 
% %      1
     

% Fit the model
ModelOutput = ienv(X, Y, u3)

% % ModelOutput = 
% % 
% %          beta: [4x2 double]
% %         Sigma: [4x4 double]
% %        Gamma1: [4x1 double]
% %        Gamma0: [4x3 double]
% %             B: [3x1 double]
% %          eta1: [2x1 double]
% %          eta2: [2x1 double]
% %        Omega1: 8.3751
% %        Omega0: [3x3 double]
% %         alpha: [4x1 double]
% %      paramNum: 16
% %             l: -1.4805e+03
% %     covMatrix: [8x8 double]
% %       asyIenv: [4x2 double]
% %         ratio: [4x2 double]
% %             n: 150


ModelOutput.beta % Print estimated regression coefficients

% % ans =
% % 
% %   -15.7396   -6.0164
% %     4.5547   -1.9483
% %   -40.8698  -12.7309
% %   -17.7939   -6.9617
   
   
ModelOutput.ratio % Check the relative efficiency of the inner envelope model versus the standard model

% % ans =
% % 
% %     1.0049    1.2694
% %     1.0020    1.0876
% %     1.0034    1.1634
% %     1.0004    1.0140
    

eig(ModelOutput.Omega1) % Check the eigenvalues of \Omega

% % ans =
% % 
% %     8.3751
    
    
eig(ModelOutput.Omega0) % Check the eigenvalues of \Omega_0

% % ans =
% % 
% %    43.5202
% %     5.5011
% %     2.1928

    
% Inference tools
B = 100;
bootse = bootstrapse(X, Y, u3, B, 'ienv') % Compute bootstrap standard errors for elements in \beta

% % bootse =
% % 
% %    13.1213    4.7810
% %     7.1204    2.6493
% %    14.5164    5.1028
% %     8.0796    2.9061


Xnew = mean(X)';
PredictOutput = prediction(ModelOutput, Xnew, 'prediction', 'ienv') % Prediction based on a new observation

% % PredictOutput = 
% % 
% %         value: [4x1 double]
% %     covMatrix: [4x4 double]
% %            SE: [4x1 double]


PredictOutput.value % Check predicted value

% % ans =
% % 
% %    58.4333
% %    30.5733
% %    37.5800
% %    11.9933


PredictOutput.SE % Check the standard errors of each elements in the predicted value

% % ans =
% % 
% %     5.1480
% %     3.3450
% %     4.3405
% %     2.0376


TestOutput = testcoefficient(ModelOutput, 'ienv'); % Test if \beta = 0

% % Test Hypothesis     Chisq Statistic    DF     P-value
% % -------------------------------------------------------
% % L * beta * R = A        4642.913        8      0.0000
% % -------------------------------------------------------


TestOutput

% % TestOutput = 
% % 
% %     chisqStatistic: 4.6429e+03
% %                 df: 8
% %             pValue: 0
% %          covMatrix: [8x8 double]


TestInput.L = ones(1, 4); % Define L
TestInput.R = [-1 2]';
TestInput.A = 0; % Define A
TestOutput = testcoefficient(ModelOutput, 'ienv', TestInput); % Test if L * \beta = A

% % Test Hypothesis     Chisq Statistic    DF     P-value
% % -------------------------------------------------------
% % L * beta * R = A          17.954        1      0.0000
% % -------------------------------------------------------


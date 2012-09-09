% This example demonstrates the usage of envmean module.

load Adopted  % Load data
Y = Adopted(:, 1 : 6); % Assign the variables


% Dimension selection
u1 = aic_envmean(Y) % Select u using AIC

% % u1 =
% % 
% %      3


u2 = bic_envmean(Y) % Select u using BIC

% % u2 =
% % 
% %      3

     
alpha = 0.01;
u3 = lrt_envmean(Y, alpha) % Select u using LRT

% % u3 =
% % 
% %      2
     

% Fit the model
ModelOutput = envmean(Y, u1)

% % ModelOutput = 
% % 
% %            mu: [6x1 double]
% %         Sigma: [6x6 double]
% %         Gamma: [6x3 double]
% %        Gamma0: [6x3 double]
% %           eta: [3x1 double]
% %         Omega: [3x3 double]
% %        Omega0: [3x3 double]
% %             l: -1.3492e+03
% %     covMatrix: [6x6 double]
% %        asyEnv: [6x1 double]
% %         ratio: [6x1 double]
% %      paramNum: 24
% %             n: 62


ModelOutput.mu % Print estimated multivariate mean

% % ans =
% % 
% %    12.3258
% %    85.9841
% %   115.5767
% %   112.1291
% %   114.4862
% %   106.4240
   
   
ModelOutput.Sigma % Check estimated error covariance matrix

% % ans =
% % 
% %     8.3278    2.9150   -4.0008    0.1057    0.3731   -0.8157
% %     2.9150  235.6587    5.8146   42.1613   59.7492   71.9394
% %    -4.0008    5.8146  179.2066   91.7228   92.9563   72.1815
% %     0.1057   42.1613   91.7228  167.1073  114.4203  110.2918
% %     0.3731   59.7492   92.9563  114.4203  184.3248  161.8752
% %    -0.8157   71.9394   72.1815  110.2918  161.8752  233.9185
    

eig(ModelOutput.Omega) % Check the eigenvalues of \Omega

% % ans =
% % 
% %     8.1096
% %   550.7479
% %   120.9889
    
    
eig(ModelOutput.Omega0) % Check the eigenvalues of \Omega_0

% % ans =
% % 
% %   220.4460
% %    69.4019
% %    38.8495

    
% Inference tools
B = 100;
bootse = bstrp_envmean(Y, u1, B) % Compute bootstrap standard errors for elements in \mu

% % bootse =
% % 
% %     1.1188
% %    11.8362
% %    13.3159
% %    18.7184
% %    25.4904
% %    26.6262


PredictOutput = predict_envmean(ModelOutput, 'prediction') % Predict a future observation

% % PredictOutput = 
% % 
% %         value: [6x1 double]
% %     covMatrix: [6x6 double]
% %            SE: [6x1 double]


PredictOutput.value % Check predicted value

% % ans =
% % 
% %    12.3258
% %    85.9841
% %   115.5767
% %   112.1291
% %   114.4862
% %   106.4240


PredictOutput.SE % Check the standard errors of each elements in the predicted value

% % ans =
% % 
% %     2.9090
% %    15.4745
% %    13.4943
% %    13.0308
% %    13.6857
% %    15.4172


TestOutput = testcoefficient_envmean(ModelOutput); % Test if \mu = 0

% %  Test Hypothesis     Chisq Statistic    DF     P-value
% % -------------------------------------------------------
% % L * mu = A             8716.691          6      0.0000
% % -------------------------------------------------------


TestOutput

% % TestOutput = 
% % 
% %     chisqStatistic: 8.7167e+03
% %                 df: 6
% %             pValue: 0
% %          covMatrix: [6x6 double]


TestInput.L = ones(1, 6); % Define L
TestInput.A = 0; % Define A
TestOutput = testcoefficient_envmean(ModelOutput, TestInput); % Test if L * \mu = A

% %  Test Hypothesis     Chisq Statistic    DF     P-value
% % -------------------------------------------------------
% % L * mu = A             6993.498          1      0.0000
% % -------------------------------------------------------



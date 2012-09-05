% This example demonstrates the usage of envmean module.

load wheatprotein.txt  % Load data
Y = wheatprotein(:, 1 : 6); % Assign the variables


% Dimension selection
u1 = aic_envmean(Y) % Select u using AIC

% % u1 =
% % 
% %      6


u2 = bic_envmean(Y) % Select u using BIC

% % u2 =
% % 
% %      5

     
alpha = 0.01;
u3 = lrt_envmean(Y, alpha) % Select u using LRT

% % u3 =
% % 
% %      5
     

% Fit the model
ModelOutput = envmean(Y, u3)

% % ModelOutput = 
% % 
% %            mu: [6x1 double]
% %         Sigma: [6x6 double]
% %         Gamma: [6x5 double]
% %        Gamma0: [6x1 double]
% %           eta: [5x1 double]
% %         Omega: [5x5 double]
% %        Omega0: 214.1646
% %             l: -881.6301
% %     covMatrix: [6x6 double]
% %        asyEnv: [6x1 double]
% %         ratio: [6x1 double]
% %            np: 26
% %             n: 50


ModelOutput.mu % Print estimated multivariate mean

% % ans =
% % 
% %   474.0794
% %   129.7194
% %   253.0441
% %   377.6816
% %   381.4752
% %    -7.2039
   
   
ModelOutput.ratio % Check estimated error covariance matrix

% % ans =
% % 
% %    1.0e+03 *
% % 
% %   Columns 1 through 4
% % 
% %     1.2430    1.0297    1.1037    1.1862
% %     1.0297    0.8668    0.9252    0.9733
% %     1.1037    0.9252    0.9893    1.0459
% %     1.1862    0.9733    1.0459    1.1463
% %     1.5078    1.2439    1.3413    1.4748
% %     0.6631    0.5555    0.5926    0.6382
% % 
% %   Columns 5 through 6
% % 
% %     1.5078    0.6631
% %     1.2439    0.5555
% %     1.3413    0.5926
% %     1.4748    0.6382
% %     2.1477    0.8132
% %     0.8132    0.3788
    

eig(ModelOutput.Omega) % Check the eigenvalues of \Omega

% % ans =
% % 
% %    1.0e+03 *
% % 
% %     6.5107
% %     0.0261
% %     0.0203
% %     0.0003
% %     0.0004
    
    
eig(ModelOutput.Omega0) % Check the eigenvalues of \Omega_0

% % ans =
% % 
% %   214.1646

    
% Inference tools
B = 100;
bootse = bstrp_envmean(Y, u3, B) % Compute bootstrap standard errors for elements in \mu

% % bootse =
% % 
% %     5.0647
% %     4.2156
% %     4.5189
% %     4.9283
% %     7.3856
% %     2.7329


PredictOutput = predict_envmean(ModelOutput, 'prediction') % Predict a future observation

% % PredictOutput = 
% % 
% %         value: [6x1 double]
% %     covMatrix: [6x6 double]
% %            SE: [6x1 double]


PredictOutput.value % Check predicted value

% % ans =
% % 
% %   474.0794
% %   129.7194
% %   253.0441
% %   377.6816
% %   381.4752
% %    -7.2039


PredictOutput.SE % Check the standard errors of each elements in the predicted value

% % ans =
% % 
% %    35.6077
% %    29.7337
% %    31.7655
% %    34.1944
% %    46.8043
% %    19.6574


TestOutput = testcoefficient_envmean(ModelOutput); % Test if \mu = 0

% %  Test Hypothesis     Chisq Statistic    DF     P-value
% % -------------------------------------------------------
% % L * mu = A          1351183.290          6      0.0000
% % -------------------------------------------------------


TestOutput

% % TestOutput = 
% % 
% %     chisqStatistic: 1.3512e+06
% %                 df: 6
% %             pValue: 0
% %          covMatrix: [6x6 double]


TestInput.L = ones(1, 6); % Define L
TestInput.A = 0; % Define A
TestOutput = testcoefficient_envmean(ModelOutput, TestInput); % Test if L * \mu = A

% %  Test Hypothesis     Chisq Statistic    DF     P-value
% % -------------------------------------------------------
% % L * mu = A             3501.358          1      0.0000
% % -------------------------------------------------------



% This example demonstrates the usage of henv module.
% The data is used as an example in Su and Cook (2012).

load waterstrider.mat  % Load data
unique(X, 'rows')  % Take a look at X.  With henv, X is group indicator which takes p different values for p groups.

% % ans =
% % 
% %     -1    -1
% %      0     1
% %      1     0


% Dimension selection
u1 = modelselectaic(X, Y, 'henv') % Select u using AIC

% % u1 =
% % 
% %      6


u2 = modelselectbic(X, Y, 'henv') % Select u using BIC

% % u2 =
% % 
% %      4

     
alpha = 0.01;
u3 = modelselectlrt(X, Y, alpha, 'henv') % Select u using LRT

% % u3 =
% % 
% %      6


% Fit the model
ModelOutput = henv(X, Y, u3)

% % ModelOutput = 
% % 
% %            mu: [8x1 double]
% %           mug: [8x3 double]
% %          Yfit: [90x8 double]
% %         Gamma: [8x6 double]
% %        Gamma0: [8x2 double]
% %          beta: [8x3 double]
% %      groupInd: [3x2 double]
% %         Sigma: [8x8x3 double]
% %           eta: [6x3 double]
% %         Omega: [6x6x3 double]
% %        Omega0: [2x2 double]
% %      paramNum: 98
% %             l: 1.0051e+03
% %     covMatrix: [32x32 double]
% %       asyHenv: [8x3 double]
% %         ratio: [8x3 double]
% %            ng: [3x1 double]


ModelOutput.mug  % Check the mean for the three groups.  Each column of the output corresponds to a group mean.

% % ans =
% % 
% %    -1.1417   -1.1267   -1.0845
% %    -1.4063   -1.4067   -1.3132
% %    -1.3314   -1.3336   -1.2152
% %    -0.3113   -0.1839   -0.1736
% %     0.4003    0.3847    0.3072
% %     0.4107    0.3753    0.3735
% %     0.3467    0.3271    0.3179
% %    -0.1954   -0.2100   -0.3488


ModelOutput.Sigma % Check the error covairance matrix.  This is a three dimensional matrix, with the ith depth containing the error covairance matrix of the ith group.

% % ModelOutput.Sigma
% % 
% % ans(:,:,1) =
% % 
% %   Columns 1 through 4
% % 
% %     0.1335    0.1302    0.1226    0.0825
% %     0.1302    0.1349    0.1251    0.0847
% %     0.1226    0.1251    0.1178    0.0790
% %     0.0825    0.0847    0.0790    0.0541
% %     0.1611    0.1637    0.1530    0.1037
% %     0.1233    0.1253    0.1171    0.0793
% %     0.1657    0.1677    0.1567    0.1061
% %     0.1358    0.1371    0.1281    0.0863
% % 
% %   Columns 5 through 8
% % 
% %     0.1611    0.1233    0.1657    0.1358
% %     0.1637    0.1253    0.1677    0.1371
% %     0.1530    0.1171    0.1567    0.1281
% %     0.1037    0.0793    0.1061    0.0863
% %     0.2047    0.1556    0.2084    0.1692
% %     0.1556    0.1193    0.1594    0.1301
% %     0.2084    0.1594    0.2135    0.1741
% %     0.1692    0.1301    0.1741    0.1439
% % 
% % 
% % ans(:,:,2) =
% % 
% %   Columns 1 through 4
% % 
% %     0.1336    0.1308    0.1222    0.0821
% %     0.1308    0.1383    0.1248    0.0827
% %     0.1222    0.1248    0.1209    0.0779
% %     0.0821    0.0827    0.0779    0.0545
% %     0.1609    0.1635    0.1524    0.1046
% %     0.1232    0.1243    0.1166    0.0800
% %     0.1658    0.1672    0.1563    0.1069
% %     0.1360    0.1364    0.1284    0.0868
% % 
% %   Columns 5 through 8
% % 
% %     0.1609    0.1232    0.1658    0.1360
% %     0.1635    0.1243    0.1672    0.1364
% %     0.1524    0.1166    0.1563    0.1284
% %     0.1046    0.0800    0.1069    0.0868
% %     0.2043    0.1559    0.2083    0.1697
% %     0.1559    0.1198    0.1596    0.1301
% %     0.2083    0.1596    0.2137    0.1740
% %     0.1697    0.1301    0.1740    0.1433
% % 
% % 
% % ans(:,:,3) =
% % 
% %   Columns 1 through 4
% % 
% %     0.1338    0.1315    0.1219    0.0820
% %     0.1315    0.1405    0.1231    0.0834
% %     0.1219    0.1231    0.1240    0.0791
% %     0.0820    0.0834    0.0791    0.0556
% %     0.1608    0.1630    0.1518    0.1029
% %     0.1233    0.1247    0.1167    0.0799
% %     0.1655    0.1663    0.1559    0.1064
% %     0.1359    0.1361    0.1279    0.0871
% % 
% %   Columns 5 through 8
% % 
% %     0.1608    0.1233    0.1655    0.1359
% %     0.1630    0.1247    0.1663    0.1361
% %     0.1518    0.1167    0.1559    0.1279
% %     0.1029    0.0799    0.1064    0.0871
% %     0.2058    0.1557    0.2088    0.1699
% %     0.1557    0.1198    0.1595    0.1300
% %     0.2088    0.1595    0.2147    0.1742
% %     0.1699    0.1300    0.1742    0.1436


ModelOutput.ratio % Check the relative efficiency of the heteroscedastic envelope model versus the standard model

% % ans =
% % 
% %     6.5439   11.2830    6.4954
% %     4.6325    5.3226    4.7242
% %     4.4456    5.0741    4.4198
% %     4.7338    6.2469    5.1937
% %     8.0377   12.5386    9.4823
% %     9.5067   11.5974   11.3444
% %    11.8632   15.6080   12.5611
% %     6.9792   11.1559   10.1002


% Inference tools
B = 100;
bootse = bootstrapse(X, Y, u3, B, 'henv') % Compute bootstrap standard errors for elements in \beta

% % bootse =
% % 
% %     0.0355    0.0459    0.0710
% %     0.0368    0.0488    0.0756
% %     0.0358    0.0405    0.0675
% %     0.0230    0.0297    0.0467
% %     0.0442    0.0569    0.0895
% %     0.0339    0.0436    0.0691
% %     0.0448    0.0583    0.0918
% %     0.0360    0.0474    0.0738


Xnew = [-1 -1]';
PredictOutput = prediction(ModelOutput, Xnew, 'prediction', 'henv') % Prediction for a new observation given the group indicator.

% % PredictOutput = 
% % 
% %         value: [8x1 double]
% %     covMatrix: [8x8 double]
% %            SE: [8x1 double]


PredictOutput.value % Check predicted value.  The predicted value is the estimated mean of the corresponding group.

% % ans =
% % 
% %    -1.1417
% %    -1.4063
% %    -1.3314
% %    -0.3113
% %     0.4003
% %     0.4107
% %     0.3467
% %    -0.1954


PredictOutput.SE % Check the standard errors of each elements in the predicted value

% % ans =
% % 
% %     0.3716
% %     0.3741
% %     0.3496
% %     0.2368
% %     0.4601
% %     0.3512
% %     0.4698
% %     0.3857


TestOutput = testcoefficient(ModelOutput, 'henv'); % Test if \beta = 0, which means that all the groups have the same mean.

% % Test Hypothesis     Chisq Statistic    DF     P-value
% % -------------------------------------------------------
% % L * beta * R = A         226.256       16      0.0000
% % -------------------------------------------------------


TestOutput

% % TestOutput = 
% % 
% %     chisqStatistic: 226.2561
% %                 df: 16
% %             pValue: 0
% %          covMatrix: [16x16 double]


TestInput.R = [1 1 -2]';
TestInput.A = zeros(8, 1); % Define A
TestOutput = testcoefficient(ModelOutput, 'henv', TestInput); % Test if \beta * R = A, this is a contrast if the sum of the means of the first two groups equal to twice the mean of the third group.

% % Test Hypothesis     Chisq Statistic    DF     P-value
% % -------------------------------------------------------
% % L * beta * R = A         144.509        8      0.0000
% % -------------------------------------------------------


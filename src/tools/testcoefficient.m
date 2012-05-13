%% testcoefficient
% 
%  This function tests the null hypothesis L * beta * R = A versus the
%  alternative hypothesis L * beta * R ~= A, where beta is estimated under
%  the model in the envelope family.

%% Syntax
%         TestOutput = testcoefficient(ModelOutput, modelType) 
%         TestOutput = testcoefficient(ModelOutput, modelType, TestInput)
% 
%% Input
% 
% *ModelOutput*: A list containing the model outputs from fitting the models.
% 
% *modelType*: A string characters indicting the model, choices can be 'env',
% 'henv', 'ienv', 'penv', 'senv' and 'xenv'.
% 
% *TestInput*: A list that specifies the null hypothesis, including L, R, and
% A.  If not provided by the user, default values will be used.
%
% * TestInput.L: The matrix multiplied to $$\beta$ on the left.  According 
% to different model, it has different size requirement.  Default value will be
% set if the user does not specify.
% 
% * TestInput.R: The matrix multiplied to $$\beta$ on the right.  According 
% to different model, it has different size requirement.  Default value will be
% set if the user does not specify.
% 
% * TestInput.A: The matrix on the right handside of the equation.  Default value will be
% set if the user does not specify.
% 
%% Output
% 
% *TestOutput*: A list containing test statistics, degrees of freedom for the
% reference chi-squared distribution, and the p-value.  At the same time, a
% table is printed out.
% 
% * TestOutput.chisqStatistic: The test statistics. A real number. 
% * TestOutput.df: The degrees of freedom of the reference chi-squared
% distribution.  A positive integer.
% * TestOutput.pValue: p-value of the test.  A real number in [0, 1].

%% Description
% This function tests for hypothesis $$H_0: L\beta R = A$, versus $$H_\alpha:
% L\beta R\neq A$.  The $$\beta$ is estimated by a model in the envelope model.  If
% the user does not specify the values for L, R and A, then the test is
% equivalent to the standard F test on if $$\beta = 0$ (for 'env', 'ienv', 
% % 'penv', 'senv' and 'xenv'), or if the group main effects are all zeros 
% (for 'henv').  The test statistics used is vec $$(L\beta R - A)$
% $$\hat{\Sigma}^{-1}$ vec $$(L\beta R - A)^{T}$, and the reference
% distribution is chi-squared distribution with degrees of freedom the same
% as the length of vec(A).

%% Example
%         load wheatprotein.txt
%         X = wheatprotein(:, 8);
%         Y = wheatprotein(:, 1:6);
%         alpha = 0.01;
%         u = lrt_env(X, Y, alpha)
%         ModelOutput = env(X, Y, u)
%         modelType = 'env';
%         TestOutout = testcoefficient(ModelOutput, modelType);
% 
%         load fiberpaper.dat
%         Y = fiberpaper(:, 1 : 4);
%         Xtemp = fiberpaper(:, 5 : 7);
%         X.X1 = Xtemp(:, 3);
%         X.X2 = Xtemp(:, 1 : 2);
%         alpha = 0.01;
%         u = lrt_penv(X, Y, alpha)
%         ModelOutput = penv(X, Y, u)
%         r = size(Y, 2);
%         p1 = size(X.X1, 2);
%         TestInput.L = rand(2, r);
%         TestInput.R = rand(p1, 1);
%         TestInput.A = zeros(2, 1);
%         TestOutout = testcoefficient_penv(ModelOutput, TestInput); 

function TestOutput = testcoefficient(ModelOutput, modelType, TestInput)

% Verify and initialize the parameters
%
if nargin < 2
    
    error('Inputs: ModelOutput and modelType should be specified!');
    
elseif nargin == 2
    
    switch(modelType)
        case 'env'
            TestOutput = testcoefficient_env(ModelOutput);
        case 'henv'
            TestOutput = testcoefficient_henv(ModelOutput);
        case 'ienv'
            TestOutput = testcoefficient_ienv(ModelOutput);
        case 'penv'
            TestOutput = testcoefficient_penv(ModelOutput);
        case 'senv'
            TestOutput = testcoefficient_senv(ModelOutput);
        case 'xenv'
            TestOutput = testcoefficient_xenv(ModelOutput);
        otherwise
            fprintf('The value specified in modelType is not supported!');
    end
    
elseif nargin == 3
    
    switch(modelType)
        case 'env'
            TestOutput = testcoefficient_env(ModelOutput, TestInput);
        case 'henv'
            TestOutput = testcoefficient_henv(ModelOutput, TestInput);
        case 'ienv'
            TestOutput = testcoefficient_ienv(ModelOutput, TestInput);
        case 'penv'
            TestOutput = testcoefficient_penv(ModelOutput, TestInput);
        case 'senv'
            TestOutput = testcoefficient_senv(ModelOutput, TestInput);
        case 'xenv'
            TestOutput = testcoefficient_xenv(ModelOutput, TestInput);
        otherwise
            fprintf('The value specified in modelType is not supported!');
    end
    
end


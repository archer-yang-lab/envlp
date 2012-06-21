%% testcoefficient_envmean
% 
%  This function tests the null hypothesis L * mu = A versus the
%  alternative hypothesis L * mu ~= A, where mu is the envelope estimator
%  of the multivariate mean.

%% Syntax
%         TestOutput = testcoefficient_envmean(ModelOutput) 
%         TestOutput = testcoefficient_envmean(ModelOutput, TestInput)
% 
%% Input
% 
% *ModelOutput*: A list containing the maximum likelihood estimators and other
% statistics inherited from envmean.
% 
% *TestInput*: A list that specifies the null hypothesis, including L and
% A.  If not provided by the user, default values will be used.
%
% * TestInput.L: The matrix multiplied to $$\mu$ on the left.  It is a d1
% by p matrix, while d1 is less than or equal to p.  Default value:
% identity matrix $$I_p$.
% * TestInput.A: The vector on the right hand side of the equation.  It is a
% d1 dimensional column vector.  Default value: d1 by d2 zero matrix.
% 
%% Output
% 
% *TestOutput*: A list containing test statistics, degrees of freedom for the
% reference chi-squared distribution, the p-value, and the covariance matrix 
% of L$$\mu$.  At the same time, a table is printed out.
% 
% * TestOutput.chisqStatistic: The test statistics. A real number. 
% * TestOutput.df: The degrees of freedom of the reference chi-squared
% distribution.  A positive integer.
% * TestOutput.pValue: p-value of the test.  A real number in [0, 1].
% * TestOutput.covMatrix: The covariance matrix of L$$\mu$. A d1
% dimensional column vector.

%% Description
% This function tests for hypothesis $$H_0: L\mu = A$, versus $$H_\alpha:
% L\mu\neq A$.  The $$\mu$ is estimated by the envelope model.  If
% the user does not specify the values for L and A, then the test is
% equivalent to the standard F test on if $$\mu = 0$.  The test statistics
% used is $$(L\mu - A)$ $$\hat{\Sigma}^{-1}$ $$(L\mu - A)^{T}$,
% and the reference distribution is chi-squared distribution with degrees of
% freedom d1. 

%% Example
%         load wheatprotein.txt
%         X = wheatprotein(:, 1 : 6);
%         alpha = 0.01;
%         u = lrt_envmean(X, alpha);
%         ModelOutput = envmean(X, u);
%         TestOutout = testcoefficient_envmean(ModelOutput);
%         p = size(X, 2);
%         TestInput.L = rand(2, p);
%         TestInput.A = zeros(2, 1);
%         TestOutout = testcoefficient_envmean(ModelOutput, TestInput);

function TestOutput = testcoefficient_envmean(ModelOutput, TestInput)

if nargin < 1
    
    error('Inputs: ModelOutput should be specified!');

elseif nargin == 1
    
    u = size(ModelOutput.Gamma, 2);
    if u == 0
        error('mu is a zero vector, no test is interesting.');
    end
    p = size(ModelOutput.mu, 1);
    TestInput.L = eye(p);
    TestInput.A = zeros(p, 1);
    Ls1 = p;
    
elseif nargin == 2
    
    u = size(ModelOutput.Gamma, 2);
    if u == 0
        error('mu is a zero vector, no test is interesting.');
    end
    p = size(ModelOutput.mu, 1);    
    if isfield(TestInput, 'L')
        [Ls1 Ls2] = size(TestInput.L);
        if Ls1 > p || Ls2 ~= p
            error('The size of L is not supported.')
        end
    else
        TestInput.L = eye(p);
        Ls1 = p;
    end

    
    if isfield(TestInput, 'A')
        [As1 As2] = size(TestInput.A);
        if As1 ~= Ls1 || As2 ~= 1
            error('The size of A should be the same as L * mu.')
        end
    else
        TestInput.A = zeros(Ls1, 1);
    end
    
end

L = TestInput.L;
A = TestInput.A;
mu = ModelOutput.mu;
covMatrix = ModelOutput.covMatrix;
n = ModelOutput.n;
Sigma = L * covMatrix * L' / n;
temp = L * mu - A;

chisqStatistic = temp' * inv(Sigma) * temp;
df = Ls1;
pValue = 1 - chi2cdf(chisqStatistic, df);

TestOutput.chisqStatistic = chisqStatistic;
TestOutput.df = df;
TestOutput.pValue = pValue;
TestOutput.covMatrix = Sigma;



fprintf('\n Test Hypothesis     Chisq Statistic    DF     P-value\n');
fprintf('-------------------------------------------------------\n');
fprintf('%s %20.3f   %8d  %10.4f\n', 'L * mu = A', chisqStatistic, df, pValue);
fprintf('-------------------------------------------------------\n');
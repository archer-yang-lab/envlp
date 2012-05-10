%% testcoefficient_env
% 
%  This function tests the null hypothesis L * beta * R = A versus the
%  alternative hypothesis L * beta * R ~= A, where beta is estimated under
%  the envelope model.

%% Syntax
%         TestOutput = testcoefficient_env(ModelOutput) 
%         TestOutput = testcoefficient_env(ModelOutput, TestInput)
% 
%% Input
% 
% *ModelOutput*: A list containing the maximum likelihood estimators and other
% statistics inherted from env.
% 
% *TestInput*: A list that specifies the null hypothesis, including L, R, and
% A.  If not provided by the user, default values will be used.
%
% * TestInput.L: The matrix multiplied to $$\beta$ on the left.  It is a d1
% by r matrix, while d1 is less than or equal to r.  Default value:
% identity matrix $$I_r$.
% * TestInput.R: The matrix multiplied to $$\beta$ on the right.  It is a p
% by d2 matrix, while d2 is less than or equal to p.  Default value:
% identity matrix $$I_p$.
% * TestInput.A: The matrix on the right handside of the equation.  It is a
% d1 by d2 matrix.  Default value: d1 by d2 zero matrix.
% 
%% Output
% 
% *TestOutput*: A list containing test statistics, degrees of freedom for the
% reference chi-squared distribution, and the p-value.  At the same time, a
% table is printed out.
% 
% * TestOutput.testStatistic: The test statistics. A real number. 
% * TestOutput.df: The degrees of freedom of the reference chi-squared
% distribution.  A positive integer.
% * TestOutput.pValue: p-value of the test.  A real number in [0, 1].

%% Description
% This function tests for hypothesis $$H_0: L\beta R = A$, versus $$H_\alpha:
% L\beta R\neq A$.  The $$\beta$ is estimated by the envelope model.  If
% the user does not specify the values for L, R and A, then the test is
% equivalent to the standard F test on if $$\beta = 0$.  The test statistics
% used is vec $$(L\beta R - A)$ $$\hat{\Sigma}^{-1}$ vec $$(L\beta R - A)^{T}$,
% and the reference distribution is chi-squared distribution with degrees of
% freedom d1 * d2. 

%% Example
%         load wheatprotein.txt
%         X = wheatprotein(:, 8);
%         Y = wheatprotein(:, 1:6);
%         alpha = 0.01;
%         u = lrt_env(X, Y, alpha)
%         ModelOutput = env(X, Y, u)
%         TestOutout = testcoefficient_env(ModelOutput);
%         r = size(Y, 2);
%         p = size(X, 2);
%         TestInput.L = rand(2, r);
%         TestInput.R = rand(p, 1);
%         TestInput.A = zeros(2, 1);
%         TestOutout = testcoefficient_senv(ModelOutput, TestInput); 


function TestOutput = testcoefficient_env(ModelOutput, TestInput)

if nargin < 1
    
    error('Inputs: ModelOutput should be specified!');

elseif nargin == 1
    
    [r p] = size(ModelOutput.beta);
    TestInput.L = eye(r);
    TestInput.R = eye(p);
    TestInput.A = zeros(r, p);
    Ls1 = r;
    Rs2 = p;
    
elseif nargin == 2
    
    [r p] = size(ModelOutput.beta);
    
    if isfield(TestInput, 'L')
        [Ls1 Ls2] = size(TestInput.L);
        if Ls1 > r || Ls2 ~= r
            error('The size of L is not supported.')
        end
    else
        TestInput.L = eye(r);
        Ls1 = r;
    end
    
    if isfield(TestInput, 'R')
        [Rs1 Rs2] = size(TestInput.R);
        if Rs1 ~= p || Rs2 > p
            error('The size of R is not supported.')
        end
    else
        TestInput.R = eye(p);
        Rs2 = p;
    end
    
    if isfield(TestInput, 'A')
        [As1 As2] = size(TestInput.A);
        if As1 ~= Ls1 || As2 ~= Rs2
            error('The size of A should be the same as L * beta * R.')
        end
    else
        TestInput.A = zeros(Ls1, Rs2);
    end
    
end

L = TestInput.L;
R = TestInput.R;
A = TestInput.A;
beta = ModelOutput.beta;
covMatrix = ModelOutput.covMatrix;
n = ModelOutput.n;

Sigma = kron(R', L) * covMatrix * kron(R, L') / n;
temp = reshape(L * beta * R - A, 1, Ls1 * Rs2);

testStatistic = temp * inv(Sigma) * temp';
df = Ls1 * Rs2;
pValue = 1 - chi2cdf(testStatistic, df);

TestOutput.testStatistic = testStatistic;
TestOutput.df = df;
TestOutput.pValue = pValue;


fprintf('\n Test Hypothesis     Test Statistics    Degrees of Freedom     P-value\n') ;
fprintf('----------------------------------------------------------------------\n') ;
fprintf('%s %13.3f   %15d  %19.4f\n', 'L * beta * R = A', testStatistic, df, pValue);
fprintf('----------------------------------------------------------------------\n') ;
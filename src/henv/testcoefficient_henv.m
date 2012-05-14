%% testcoefficient_henv
% 
%  This function tests the null hypothesis L * beta * R = A versus the
%  alternative hypothesis L * beta * R ~= A, where beta is estimated under
%  the heteroscedastic envelope model.

%% Syntax
%         TestOutput = testcoefficient_henv(ModelOutput) 
%         TestOutput = testcoefficient_henv(ModelOutput, TestInput)
% 
%% Input
% 
% *ModelOutput*: A list containing the maximum likelihood estimators and other
% statistics inherited from henv.
% 
% *TestInput*: A list that specifies the null hypothesis, including L, R, and
% A.  If not provided by the user, default values will be used.
%
% * TestInput.L: The matrix multiplied to $$\beta$ on the left.  It is a d1
% by r matrix, while d1 is less than or equal to r - 1.  Default value:
% identity matrix $$I_{r}$.  
% * TestInput.R: The matrix multiplied to $$\beta$ on the right.  It is a p
% by d2 matrix, while d2 is less than or equal to p.  Default value:
% identity matrix $$(I_{p-1}, 0_{(p-1)\times 1})^{T}$.  This is because the
% columns of $$\beta$ sum to 0.  Then we cannot use $$I_p$ as default.
% * TestInput.A: The matrix on the right hand side of the equation.  It is a
% d1 by d2 matrix.  Default value: d1 by d2 zero matrix.
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
% L\beta R\neq A$.  The $$\beta$ is estimated by the heteroscedastic envelope model.  If
% the user does not specify the values for L, R and A, then the test is
% equivalent to the standard F test on if all the main group effects are 0.  The test statistics
% used is vec $$(L\beta R - A)$ $$\hat{\Sigma}^{-1}$ vec $$(L\beta R - A)^{T}$,
% and the reference distribution is chi-squared distribution with degrees of
% freedom d1 * d2. 

%% Example
%         load waterstrider.mat
%         u = lrt_henv(X, Y, 0.01);
%         ModelOutput = henv(X, Y, u);
%         TestOutout = testcoefficient_henv(ModelOutput);
%         r = size(Y, 2);
%         p = size(ModelOutput.beta, 2);
%         TestInput.L = rand(2, r);
%         TestInput.R = rand(p, 1);
%         TestInput.A = zeros(2, 1);
%         TestOutout = testcoefficient_henv(ModelOutput, TestInput); 

function TestOutput = testcoefficient_henv(ModelOutput, TestInput)

if nargin < 1
    
    error('Inputs: ModelOutput should be specified!');

elseif nargin == 1
    
    [r p] = size(ModelOutput.beta);
    TestInput.L = eye(r);
    TestInput.R = [eye(p - 1) zeros(p - 1, 1)]';
    TestInput.A = zeros(r, p - 1);
    Ls1 = r;
    Rs2 = p - 1;
    
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
        if Rs1 ~= p || Rs2 > p - 1
            error('The size of R is not supported.')
        end
    else
        TestInput.R = [eye(p - 1) zeros(p - 1, 1)]';
        Rs2 = p - 1;
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
covMatrix = ModelOutput.covMatrix(r + 1 : r * (p + 1), r + 1 : r * (p + 1));
ng = ModelOutput.ng;
multiplier = diag(kron(1./sqrt(ng), ones(r, 1)));

Sigma = kron(R', L) * multiplier * covMatrix * multiplier * kron(R, L');
temp = reshape(L * beta * R - A, 1, Ls1 * Rs2);

chisqStatistic = temp * inv(Sigma) * temp';
df = Ls1 * Rs2;
pValue = 1 - chi2cdf(chisqStatistic, df);

TestOutput.chisqStatistic = chisqStatistic;
TestOutput.df = df;
TestOutput.pValue = pValue;


fprintf('\nTest Hypothesis     Chisq Statistic    DF     P-value\n') ;
fprintf('----------------------------------------------------------------------\n') ;
fprintf('%s %15.3f   %6d  %10.4f\n', 'L * beta * R = A', chisqStatistic, df, pValue);
fprintf('----------------------------------------------------------------------\n') ;
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
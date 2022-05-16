function [result, FitnessEvaluations] = Fitness(varargin)
global ContexVector Benchmark
narginchk(4,6);
if nargin == 4
    X = varargin{1};
    Variables = varargin{2};
    ProblemType = varargin{3};
    FitnessEvaluations = varargin{4};
elseif nargin == 6
    X = varargin{1};
    Variables = varargin{2};
    ProblemType = varargin{3};
    MapIndex = varargin{4};
    DecompContext = varargin{5};
    FitnessEvaluations = varargin{6};
end
tmp = NaN(size(ContexVector));
[SolutionNumber,~] = size(X);
Y = repmat(tmp,SolutionNumber,1);
if strcmp(ProblemType, 'independent parts')
    for ii=1 : SolutionNumber
        Y(ii,Variables) = X(ii,:);
    end
    [result, FitnessEvaluations] = benchmark_func_ind(Y, FitnessEvaluations);
elseif strcmp(ProblemType, 'decomposition parts')
    for ii=1 : SolutionNumber
        solution = DecompContext.DecompContextVector;
        x = X(ii,:);
        [~, dimensions] = size(x); 
        for jj = 1 : dimensions
            solution(:, MapIndex(:, jj)) = x(:, jj);
        end
        Y(ii,Variables) = solution;
    end
    [result, FitnessEvaluations] = benchmark_func_decomp(Y, FitnessEvaluations);
elseif strcmp(ProblemType, 'combination parts')
    [~, CombinedProblemNum] = size(Variables); 
    Variables_ = [];
    for i = 1 : CombinedProblemNum
        Variables_ = [Variables_, Variables{i}];
    end
    % sort the variable index according to Benchmark.permutationMap
%     variables__ = [];
%     for i = 1 : length(Benchmark.PermutationMap)
%         for j = 1 : length(Variables_)
%             if Benchmark.PermutationMap(:, i) == Variables_(:, j)
%                  variables__ = [variables__, Variables_(:, j)];
%                 break;
%             end
%         end
%     end
    for ii=1 : SolutionNumber
        Y(ii,Variables_) = X(ii,:);
    end
    [result, FitnessEvaluations] = benchmark_func_comb(Y, CombinedProblemNum, FitnessEvaluations);
end 
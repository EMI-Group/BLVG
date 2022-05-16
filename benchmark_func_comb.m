function [result, FitnessEvaluations] = benchmark_func_comb(X, CombinedProblemNum, FitnessEvaluations)
global Benchmark MPB evals ChangeFrequency
[SolutionNumber,~] = size(X);
result = NaN(SolutionNumber,CombinedProblemNum);
for i=1 : SolutionNumber
    x = X(i,:)';
    ldim = 1;
    count = 1;
    for ii=1:Benchmark.MPBnumber
        solution = x(Benchmark.PermutationMap(ldim:ldim+MPB{ii}.Dimension-1));
        if ~isnan(solution)
            f = Basefunction(solution,ii);
            result(i, count) = f;
            count = count + 1;
        end
        ldim = ldim + MPB{ii}.Dimension;
    end
    evals = evals + 1;
    FitnessEvaluations = FitnessEvaluations + 1;
%     disp(evals);
    if ~rem(evals,ChangeFrequency)
        EnvironmentalChange;
    end
end
end
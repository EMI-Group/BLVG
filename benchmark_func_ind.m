function [result, FitnessEvaluations] = benchmark_func_ind(X, FitnessEvaluations)
global Benchmark MPB evals ChangeFrequency
[SolutionNumber,~] = size(X);
result = NaN(SolutionNumber,1);
for i=1 : SolutionNumber
    x = X(i,:)';
    fit = 0;
    ldim = 1;
    for ii=1:Benchmark.MPBnumber
        solution = x(Benchmark.PermutationMap(ldim:ldim+MPB{ii}.Dimension-1));
        f=0;
        if ~isnan(solution)
            f = Basefunction(solution,ii);
        end
        ldim = ldim + MPB{ii}.Dimension;
        fit = fit + f;
    end
    result(i) = fit;
    evals = evals + 1;
    FitnessEvaluations = FitnessEvaluations + 1;
%     disp(evals);
    if ~rem(evals,ChangeFrequency)
        EnvironmentalChange;
    end
end
end
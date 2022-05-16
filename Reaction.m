global PSO
[~, mPSONumber] = size(PSO);
for i = 1 : mPSONumber
    tmpEvals = evals;
    if strcmp(PSO{i}.ProblemType, 'independent parts')
        PSO{i} = Reaction_independent(PSO{i});
    elseif strcmp(PSO{i}.ProblemType, 'decomposition parts')
        PSO{i} = Reaction_decomposition(PSO{i});
    elseif strcmp(PSO{i}.ProblemType, 'combination parts')
        PSO{i} = Reaction_combination(PSO{i});
    end
end
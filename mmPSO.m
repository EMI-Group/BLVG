global PSO ChangeFlag evals;
[~, mPSONumber] = size(PSO);
for i = 1 : mPSONumber
    if ChangeFlag == 0
        tmpEvals = evals;
        if PSO{i}.FitnessEvaluations < PSO{i}.MaxFitnessEvaluations
            if strcmp(PSO{i}.ProblemType, 'independent parts')
                PSO{i} = mPSO_independent(PSO{i});
%                 disp([evals - tmpEvals]);
            elseif strcmp(PSO{i}.ProblemType, 'decomposition parts')
                PSO{i} = mPSO_decomposition(PSO{i});
%                 disp([evals - tmpEvals]);
            elseif strcmp(PSO{i}.ProblemType, 'combination parts')
                PSO{i} = mPSO_combination(PSO{i});
%                 disp([evals - tmpEvals]);
            end
        end
    end
end
if ChangeFlag ~= 0
    return;
end
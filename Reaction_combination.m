function PSO = Reaction_combination(PSO)
PSO.FitnessEvaluations = 0;
%% Updating Shift Severity
for i = 1 : PSO.SwarmNumber
    if i ~= PSO.FreeSwarmID
        if sum(isnan(PSO.particle(i).Gbest_past_environment))==0
            ldim = 1;
            for j = 1 : PSO.CombinedProblemNum
                PSO.particle(i).Shifts(j) = pdist2(PSO.particle(i).Gbest_past_environment(:, ldim:ldim + length(PSO.Variable{j}) - 1), PSO.particle(i).GbestPosition(:, ldim:ldim + length(PSO.Variable{j}) - 1));
                ldim = ldim + length(PSO.Variable{j}) - 1;
            end
            PSO.particle(i).ShiftSeverity(j) = 1;
        else
            for j = 1 : PSO.CombinedProblemNum
                PSO.particle(i).ShiftSeverity(j) = 1;
            end
        end
    end
end
%% Introducing diversity (all trackers except free swarm)
for i = 1 : PSO.SwarmNumber
    if i ~= PSO.FreeSwarmID
        ldim = 1;
        for j = 1 : PSO.CombinedProblemNum
            PSO.particle(i).X(:, ldim:ldim + length(PSO.Variable{j}) - 1) = repmat(PSO.particle(i).GbestPosition(:, ldim:ldim + length(PSO.Variable{j}) - 1), PSO.PopulationSize, 1)...
                + rands(PSO.PopulationSize, length(PSO.Variable{j})) * PSO.particle(i).ShiftSeverity(j);
            ldim = ldim + length(PSO.Variable{j}) - 1;
        end
        PSO.particle(i).X(1,:) = PSO.particle(i).GbestPosition;
    end
end
%% Updating memory for all swarms
for i=1 : PSO.SwarmNumber
    [PSO.particle(i).FitnessValue, PSO.FitnessEvaluations] = Fitness(PSO.particle(i).X, PSO.Variable, PSO.ProblemType, PSO.FitnessEvaluations);
    PSO.particle(i).PbestPosition = PSO.particle(i).X;
    PSO.particle(i).PbestValue = PSO.particle(i).FitnessValue;
    PSO.particle(i).Gbest_past_environment = PSO.particle(i).GbestPosition;
    [PSO.particle(i).GbestValue, GbestID] = max(PSO.particle(i).PbestValue);
    ldim = 1;
    for j = 1 : PSO.CombinedProblemNum
        PSO.particle(i).GbestPosition(:, ldim:ldim + length(PSO.Variable{j}) - 1) = PSO.particle(i).PbestPosition(GbestID(:,j),ldim:ldim + length(PSO.Variable{j}) - 1);
        ldim = ldim + length(PSO.Variable{j});
    end
end
end


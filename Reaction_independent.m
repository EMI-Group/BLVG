function PSO = Reaction_independent(PSO)
PSO.FitnessEvaluations = 0;
%% Updating Shift Severity
for i = 1 : PSO.SwarmNumber
    if i ~= PSO.FreeSwarmID
        if sum(isnan(PSO.particle(i).Gbest_past_environment))==0
            PSO.particle(i).Shifts = [PSO.particle(i).Shifts, pdist2(PSO.particle(i).Gbest_past_environment,PSO.particle(i).GbestPosition)];
            PSO.particle(i).ShiftSeverity = mean(PSO.particle(i).Shifts);
        else
            PSO.particle(i).ShiftSeverity = 1;
        end
    end
end
%% Introduce diversity
for i = 1 : PSO.SwarmNumber
    if i ~= PSO.FreeSwarmID
        PSO.particle(i).X = repmat(PSO.particle(i).GbestPosition,PSO.PopulationSize,1) + (rands(PSO.PopulationSize,PSO.Dimension) * PSO.particle(i).ShiftSeverity);
        PSO.particle(i).X(1,:) = PSO.particle(i).GbestPosition;
    end
end
%% Updating memory
for i = 1 : PSO.SwarmNumber
    PSO.particle(i).Gbest_past_environment = PSO.particle(i).GbestPosition;
    [PSO.particle(i).FitnessValue, PSO.FitnessEvaluations] = Fitness(PSO.particle(i).X, PSO.Variable, PSO.ProblemType, PSO.FitnessEvaluations);
    PSO.particle(i).PbestValue = PSO.particle(i).FitnessValue;
    PSO.particle(i).PbestPosition = PSO.particle(i).X;
    [PSO.particle(i).GbestValue, GbestID] = max(PSO.particle(i).PbestValue);
    PSO.particle(i).GbestPosition = PSO.particle(i).PbestPosition(GbestID,:);
    PSO.particle(i).DeactivationFlag = 0;   % Activation
end
end
function Pso = mPSO_combination(Pso)
global ChangeFlag;
TmpSwarmNum = Pso.SwarmNumber;
%% Sub-swarm movement
for i = 1 : Pso.SwarmNumber
    if ChangeFlag == 0 && (Pso.particle(i).DeactivationFlag == 0)  % Update active trackers and free swarm
        if Pso.FitnessEvaluations >= Pso.MaxFitnessEvaluations
            return;
        end 
        Pso.particle(i).Velocity =  Pso.x * (Pso.particle(i).Velocity + (Pso.c1 * rand(Pso.PopulationSize,Pso.Dimension).*(Pso.particle(i).PbestPosition-Pso.particle(i).X)) + (Pso.c2 * rand(Pso.PopulationSize,Pso.Dimension).*(repmat(Pso.particle(i).GbestPosition,Pso.PopulationSize,1)-Pso.particle(i).X)));
        Pso.particle(i).X = Pso.particle(i).X + Pso.particle(i).Velocity;
        for j = 1 : Pso.PopulationSize
            for k = 1 : Pso.Dimension
                if Pso.particle(i).X(j,k) > Pso.MaxCoordinate
                    Pso.particle(i).X(j,k) = Pso.MaxCoordinate;
                    Pso.particle(i).Velocity(j,k) = 0;
                elseif Pso.particle(i).X(j,k) < Pso.MinCoordinate
                    Pso.particle(i).X(j,k) = Pso.MinCoordinate;
                    Pso.particle(i).Velocity(j,k) = 0;
                end
            end
        end 
        [Pso.particle(i).FitnessValue, Pso.FitnessEvaluations] = Fitness(Pso.particle(i).X, Pso.Variable, Pso.ProblemType, Pso.FitnessEvaluations);
        tmp = Pso.particle(i).FitnessValue > Pso.particle(i).PbestValue;
        Pso.particle(i).PbestValue(tmp)  = Pso.particle(i).FitnessValue(tmp);
        ldim = 1;
        for j = 1 : Pso.CombinedProblemNum
            for k = 1 : Pso.PopulationSize
                if(tmp(k,j) == 1)
                    Pso.particle(i).PbestPosition(k, ldim:ldim+length(Pso.Variable{j})-1) = Pso.particle(i).X(k, ldim:ldim+length(Pso.Variable{j})-1);
                end
            end
            ldim = ldim + length(Pso.Variable{j});
        end
        [Pso.particle(i).GbestValue,GbestID] = max(Pso.particle(i).PbestValue);
        ldim = 1;
        for j = 1 : Pso.CombinedProblemNum
            Pso.particle(i).GbestPosition(:, ldim:ldim+length(Pso.Variable{j})-1) = ...
                Pso.particle(i).PbestPosition(GbestID(:,j),ldim:ldim+length(Pso.Variable{j})-1);
            ldim = ldim + length(Pso.Variable{j});
        end
    end
end
%% Exclusion
Exclusion = zeros(Pso.SwarmNumber, Pso.CombinedProblemNum);
ldim = 1;
for i = 1 : Pso.CombinedProblemNum
    for j = 1 : Pso.SwarmNumber - 1
        for k = j+1 : Pso.SwarmNumber
            if norm((Pso.particle(j).GbestPosition(:, ldim:ldim+length(Pso.Variable{i})-1) - Pso.particle(k).GbestPosition(:, ldim:ldim+length(Pso.Variable{i})-1)), inf) < Pso.ExclusionLimit(:,i)
                if Pso.particle(j).GbestValue(:, i) < Pso.particle(k).GbestValue(:, i)
                    Exclusion(j, i) = 1;
                else
                    Exclusion(k, i) = 1;
                end
            end
        end
    end
    ldim = ldim + length(Pso.Variable{i});
end
% Delete swarms that are totally in the exclusion areas of others
count = 0;
DeleteSwarm = [];
IntegrationFlag = 0;
for i = 1 : Pso.SwarmNumber
    if sum(Exclusion(i,:)) == Pso.CombinedProblemNum
        Pso.particle(i - count) = [];
        count = count + 1;
        Pso.SwarmNumber = Pso.SwarmNumber-1;
        Pso.FreeSwarmID = Pso.FreeSwarmID-1;
        DeleteSwarm = [DeleteSwarm; i];
        IntegrationFlag = 1;
    end
end
Exclusion(DeleteSwarm, :) = [];
% Integrate subcomponents according to max peaks number & delete redundent swarms 
[MaxPeaksNumber, MaxID] = min(sum(Exclusion));
MaxPeaksNumber = Pso.SwarmNumber - MaxPeaksNumber;
DeleteSwarm = find(Exclusion(:, MaxID));
[DeleteSwarmNumber, ~] = size(DeleteSwarm);
ldim = 1;
for i = 1 : Pso.CombinedProblemNum
    if i ~= MaxID
        for j = 1 : DeleteSwarmNumber
            ReplaceFlag = 0;
            if Exclusion(DeleteSwarm(j,:), i) == 0 && ReplaceFlag == 0  %Save
                for k = 1 : Pso.SwarmNumber
                    ind = find(DeleteSwarm == k);
                    if size(ind,1) == 0 && Exclusion(k, i) == 1
                        Exclusion(k, i) = 0;
                        %Replace Component Exclusion(k, i) by Exclusion(DeleteSwarm(j,:), i)
                        Pso.particle(k).X(:, ldim:ldim+length(Pso.Variable{i})-1) = Pso.particle(DeleteSwarm(j,:)).X(:, ldim:ldim+length(Pso.Variable{i})-1);
                        Pso.particle(k).Velocity(:, ldim:ldim+length(Pso.Variable{i})-1) = Pso.particle(DeleteSwarm(j,:)).Velocity(:, ldim:ldim+length(Pso.Variable{i})-1);
                        Pso.particle(k).FitnessValue(:, i) = Pso.particle(DeleteSwarm(j,:)).FitnessValue(:, i);
                        Pso.particle(k).PbestPosition(:, ldim:ldim+length(Pso.Variable{i})-1) = Pso.particle(DeleteSwarm(j,:)).PbestPosition(:, ldim:ldim+length(Pso.Variable{i})-1);
                        Pso.particle(k).GbestPosition(:, ldim:ldim+length(Pso.Variable{i})-1) = Pso.particle(DeleteSwarm(j,:)).GbestPosition(:, ldim:ldim+length(Pso.Variable{i})-1);
                        Pso.particle(k).PbestValue(:, i) = Pso.particle(DeleteSwarm(j,:)).PbestValue(:, i);
                        Pso.particle(k).GbestValue(:, i) = Pso.particle(DeleteSwarm(j,:)).GbestValue(:, i);
                        ReplaceFlag = 1;
                        break;
                    end
                end
            end
        end
    end 
    ldim = ldim + length(Pso.Variable{i});
end
count = 0;
for i = 1 : DeleteSwarmNumber
    Pso.particle(DeleteSwarm(i,:)-count) = [];
    Pso.SwarmNumber = Pso.SwarmNumber-1;
    Pso.FreeSwarmID = Pso.FreeSwarmID-1;
    IntegrationFlag = 1;
    count = count + 1;
end
if Pso.FitnessEvaluations >= Pso.MaxFitnessEvaluations
    return;
end
% if IntegrationFlag == 1
    Pso.SwarmNumber = Pso.SwarmNumber+1;
    Pso.FreeSwarmID = Pso.SwarmNumber;
    [Pso.particle(Pso.FreeSwarmID), Pso.FitnessEvaluations] = InitializingPSO(Pso);
% end
%% FreeSwarm Convergence Detection
if Pso.FitnessEvaluations >= Pso.MaxFitnessEvaluations
    return;
end
CovergenceFlag = 1;
for i=1 : Pso.PopulationSize
    for j=1 : Pso.PopulationSize
        if CovergenceFlag==1 && i~=j
            ldim = 1;
            for k = 1 : Pso.CombinedProblemNum
                if norm(Pso.particle(Pso.FreeSwarmID).X(i, ldim:ldim+length(Pso.Variable{k})-1)...
                        - Pso.particle(Pso.FreeSwarmID).X(j, ldim:ldim+length(Pso.Variable{k})-1), inf) > Pso.ConvergenceLimit(:,k)
                    CovergenceFlag = 0;
                    break;
                end
                ldim = ldim + length(Pso.Variable{k});
            end 
        end
    end
end
if CovergenceFlag==1
    Pso.SwarmNumber = Pso.SwarmNumber+1;
    Pso.FreeSwarmID = Pso.SwarmNumber;
    [Pso.particle(Pso.FreeSwarmID), Pso.FitnessEvaluations] = InitializingPSO(Pso);
end
%% Deactivation of converged trackers
for i = 1 : Pso.SwarmNumber
    if i ~= Pso.FreeSwarmID
        Pso.particle(i).DeactivationFlag = 1;
        for j = 1 : Pso.PopulationSize
            for k = j + 1 : Pso.PopulationSize
                if Pso.particle(i).DeactivationFlag ~= 0
                    ldim = 1;
                    for r = 1 : Pso.CombinedProblemNum
                        if (norm(Pso.particle(i).X(j, ldim:ldim+length(Pso.Variable{r})-1) - Pso.particle(i).X(k, ldim:ldim+length(Pso.Variable{r})-1), inf) > Pso.DeactivationLimit)
                            Pso.particle(i).DeactivationFlag = 0;
                            break;
                        end
                        ldim = ldim + length(Pso.Variable{r});
                    end
                end
            end
        end
    end
end
%% Updating Thresholds
% if TmpSwarmNum ~= Pso.SwarmNumber
%     for i = 1 : Pso.CombinedProblemNum
%         Pso.ExclusionLimit(i) = 0.5 * ((Pso.MaxCoordinate-Pso.MinCoordinate) / ((Pso.SwarmNumber) ^ (1 / Pso.Dimension)));
% %         if Pso.ExclusionLimit(i) > Pso.MinExclusionLimit(i)
% %             Pso.ExclusionLimit(i) = Pso.MinExclusionLimit(i);
% %         end
%     end
% end

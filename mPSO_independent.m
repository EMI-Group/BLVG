function Pso = mPSO_independent(Pso)
global ChangeFlag
TmpSwarmNum = Pso.SwarmNumber;
%%Sub-swarm movement
for ii=1 : Pso.SwarmNumber
    if ChangeFlag==0 && (Pso.particle(ii).DeactivationFlag == 0)
        if Pso.FitnessEvaluations >= Pso.MaxFitnessEvaluations
            return;
        end
        Pso.particle(ii).Velocity = Pso.x * (Pso.particle(ii).Velocity + (Pso.c1 * rand(Pso.PopulationSize , Pso.Dimension).*(Pso.particle(ii).PbestPosition - Pso.particle(ii).X)) + (Pso.c2*rand(Pso.PopulationSize , Pso.Dimension).*(repmat(Pso.particle(ii).GbestPosition,Pso.PopulationSize,1) - Pso.particle(ii).X)));
        Pso.particle(ii).X = Pso.particle(ii).X + Pso.particle(ii).Velocity;
        for jj=1 : Pso.PopulationSize
            for kk=1 : Pso.Dimension
                if Pso.particle(ii).X(jj,kk) > Pso.MaxCoordinate
                    Pso.particle(ii).X(jj,kk) = Pso.MaxCoordinate;
                    Pso.particle(ii).Velocity(jj,kk) = 0;
                elseif Pso.particle(ii).X(jj,kk) < Pso.MinCoordinate
                    Pso.particle(ii).X(jj,kk) = Pso.MinCoordinate;
                    Pso.particle(ii).Velocity(jj,kk) = 0;
                end
            end
        end 
        [Pso.particle(ii).FitnessValue, Pso.FitnessEvaluations]  = Fitness(Pso.particle(ii).X, Pso.Variable, Pso.ProblemType, Pso.FitnessEvaluations);
        for jj=1 : Pso.PopulationSize
            if Pso.particle(ii).FitnessValue(jj) > Pso.particle(ii).PbestValue(jj)
                Pso.particle(ii).PbestValue(jj) = Pso.particle(ii).FitnessValue(jj);
                Pso.particle(ii).PbestPosition(jj,:) = Pso.particle(ii).X(jj,:);
            end
        end
        [BestPbestValue,BestPbestID] = max(Pso.particle(ii).PbestValue);
        if BestPbestValue > Pso.particle(ii).GbestValue
            Pso.particle(ii).GbestValue = BestPbestValue;
            Pso.particle(ii).GbestPosition = Pso.particle(ii).PbestPosition(BestPbestID,:);
        end
    end
end
%% Exclusion
ExclusionFlag = 0;
for ii=1 : Pso.SwarmNumber
    for jj=ii+1 : Pso.SwarmNumber
        if ExclusionFlag==0 && pdist2(Pso.particle(ii).GbestPosition,Pso.particle(jj).GbestPosition)<Pso.ExclusionLimit
            if Pso.particle(ii).GbestValue<Pso.particle(jj).GbestValue
                if Pso.FreeSwarmID~=ii
                   Pso.particle(ii) = [];
                   Pso.SwarmNumber = Pso.SwarmNumber-1;
                   Pso.FreeSwarmID = Pso.FreeSwarmID-1;
                else
                   if Pso.FitnessEvaluations >= Pso.MaxFitnessEvaluations
                       return;
                   end
                   [Pso.particle(ii), Pso.FitnessEvaluations] = InitializingPSO(Pso); 
                end
            else
                if Pso.FreeSwarmID~=jj
                   Pso.particle(jj) = [];
                   Pso.SwarmNumber = Pso.SwarmNumber-1;
                   Pso.FreeSwarmID = Pso.FreeSwarmID-1;
                else
                   if Pso.FitnessEvaluations >= Pso.MaxFitnessEvaluations
                       return;
                   end
                   [Pso.particle(jj), Pso.FitnessEvaluations] = InitializingPSO(Pso); 
                end
            end
            ExclusionFlag = 1;
        end
    end
end
%% FreeSwarm Convergence
if Pso.FitnessEvaluations >= Pso.MaxFitnessEvaluations
   return;
end
CovergenceFlag = 1;
for ii=1 : Pso.PopulationSize
    for jj=1 : Pso.PopulationSize
        if CovergenceFlag==1 && ii~=jj
            if pdist2(Pso.particle(Pso.FreeSwarmID).X(ii,:),Pso.particle(Pso.FreeSwarmID).X(jj,:))>Pso.ConvergenceLimit
                CovergenceFlag = 0;
            end
        end
    end
end
if CovergenceFlag == 1
   Pso.SwarmNumber = Pso.SwarmNumber+1; 
   Pso.FreeSwarmID = Pso.SwarmNumber;
   Pso.particle(Pso.FreeSwarmID) = InitializingPSO(Pso);
end
%% Deactivation of converged trackers
for ii = 1 : Pso.SwarmNumber - 1
    Pso.particle(ii).DeactivationFlag = 1;
    for jj = 1 : Pso.PopulationSize
        for kk = jj + 1 : Pso.PopulationSize
            distPair = norm (Pso.particle(ii).X(jj, :) - Pso.particle(ii).X(kk, :), inf);
            if (distPair > Pso.DeactivationLimit)
                Pso.particle(ii).DeactivationFlag = 0;
                break;
            end
        end
    end
end
%% Updating Thresholds
if TmpSwarmNum ~= Pso.SwarmNumber
    Pso.ExclusionLimit = 0.5 * ((Pso.MaxCoordinate-Pso.MinCoordinate) / ((Pso.SwarmNumber) ^ (1 / Pso.Dimension)));
end
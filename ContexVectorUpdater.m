global PSO ContexVector
for i = 1 : length(PSO)
   if strcmp(PSO{i}.ProblemType, 'independent parts')
       MaxValue = -inf; 
       for j = 1 : PSO{i}.SwarmNumber
           if PSO{i}.particle(j).GbestValue > MaxValue
            MaxValue = PSO{i}.particle(j).GbestValue;
            MaxID = j;
           end
       end
       ContexVector(PSO{i}.Variable) = PSO{i}.particle(MaxID).GbestPosition;
   elseif strcmp(PSO{i}.ProblemType, 'decomposition parts')
       for j = 1 : PSO{i}.K
           tmp = PSO{i}.SearchSwarms{j}.particle.FitnessValue > PSO{i}.SearchSwarms{j}.particle.PbestValue;
           PSO{i}.SearchSwarms{j}.particle.PbestValue(tmp)  = PSO{i}.SearchSwarms{j}.particle.FitnessValue(tmp);
           PSO{i}.SearchSwarms{j}.particle.PbestPosition(tmp, :) = PSO{i}.SearchSwarms{j}.particle.X(tmp , :);
           [PSO{i}.SearchSwarms{j}.particle.GbestValue, PSO{i}.SearchSwarms{j}.particle.GbestPositionIndex] = max(PSO{i}.SearchSwarms{j}.particle.PbestValue);
           PSO{i}.SearchSwarms{j}.particle.GbestPosition = PSO{i}.SearchSwarms{j}.particle.PbestPosition(PSO{i}.SearchSwarms{j}.particle.GbestPositionIndex,:); 
           % Save the best results before environmental change
           if PSO{i}.SearchSwarms{j}.particle.GbestValue > PSO{i}.DecompContext.DecompContextVectorFitness
               PSO{i}.DecompContext.DecompContextVector(:, PSO{i}.SearchSwarms{j}.MapIndex) = PSO{i}.SearchSwarms{j}.particle.GbestPosition;
               DecompContextVector_ = NaN(1, size(ContexVector, 2));
               DecompContextVector_(:, PSO{i}.Variable) =  PSO{i}.DecompContext.DecompContextVector;
               [PSO{i}.DecompContext.DecompContextVectorFitness, PSO{i}.FitnessEvaluations] = benchmark_func_decomp(DecompContextVector_, PSO{i}.FitnessEvaluations); 
           end
       end
       MaxValue = PSO{i}.DecompContext.DecompContextVectorFitness;
       MaxPosition = PSO{i}.DecompContext.DecompContextVector;
       for ii = 1 : PSO{i}.TrackerSwarms.SwarmNumber
           if PSO{i}.TrackerSwarms.particle(ii).GbestValue > MaxValue
               MaxValue = PSO{i}.TrackerSwarms.particle(ii).GbestValue;
               MaxPosition = PSO{i}.TrackerSwarms.particle(ii).GbestPosition;
           end
       end
       ContexVector(PSO{i}.Variable) = MaxPosition;
   elseif strcmp(PSO{i}.ProblemType, 'combination parts')
        MaxValue = -inf;
        ldim = 1;
        for kk = 1 : size(PSO{i}.particle(1).GbestValue, 2)
            MaxValue = -inf;
            for jj=1 : PSO{i}.SwarmNumber
                if PSO{i}.particle(jj).GbestValue(:, kk) > MaxValue
                    MaxValue = PSO{i}.particle(jj).GbestValue(:, kk);
                    MaxSwarmNumberID = jj;
                end
            end
            ContexVector(:, PSO{i}.Variable{kk}) = PSO{i}.particle(MaxSwarmNumberID).GbestPosition(:, ldim:ldim+length(PSO{i}.Variable{kk})-1);
            ldim = ldim + length(PSO{i}.Variable{kk});
        end 
   end 
end



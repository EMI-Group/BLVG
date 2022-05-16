function [delta, lambda, evaluations] = ism()  
%%%% Implement Algorithm2 ISM in DG2
global Benchmark;
Dimension  = Benchmark.Dimension;
FEs = 0;
Median = (Benchmark.UB + Benchmark.LB)/2;%m 

f_archive    = NaN(Dimension, Dimension);
fhat_archive = NaN(Dimension, 1);
delta1 = NaN(Dimension, Dimension);
delta2 = NaN(Dimension, Dimension);
lambda = NaN(Dimension, Dimension);

Point1 = Benchmark.LB;
results = benchmark_func(Point1);
Point1_FitnessValue = sum(results);%f_base
FEs = FEs + 1;

counter = 0;
% prev = 0;
% prog = 0;

for i=1:Dimension-1%Line 6 in Algorithm2
    if(~isnan(fhat_archive(i)))
        Point2_FitnessValue = fhat_archive(i);
    else
        Point2 = Point1;
        Point2(i) = Median(i);
        results =  benchmark_func(Point2);
        Point2_FitnessValue = sum(results);
        FEs = FEs + 1;
        fhat_archive(i) = Point2_FitnessValue;
    end
    for j=i+1:Dimension%Line 10 in Algorithm2
        if ~mod(FEs,10000)
            x = floor(counter/(Dimension*(Dimension-1))*2*100);
            clc;
            disp(['Percentage  ' , num2str(x),'%']);
        end
        counter = counter + 1;
        
        if(~isnan(fhat_archive(j)))
            Point3_FitnessValue = fhat_archive(j);
        else
            Point3 = Point1;
            Point3(j) = Median(j);
            results = benchmark_func(Point3);
            Point3_FitnessValue = sum(results);
            FEs = FEs + 1;
            fhat_archive(j) = Point3_FitnessValue;
        end
        Point4 = Point1;
        Point4(i) = Median(i);
        Point4(j) = Median(j);
        results = benchmark_func(Point4);
        Point4_FitnessValue = sum(results);
        FEs = FEs + 1;
        f_archive(i, j) = Point4_FitnessValue;
        f_archive(j, i) = Point4_FitnessValue;
        
        d1 = Point2_FitnessValue - Point1_FitnessValue;
        d2 = Point4_FitnessValue - Point3_FitnessValue;
        
        delta1(i, j) = d1;
        delta2(i, j) = d2;
        lambda(i, j) = abs(d1 - d2);
        
    end
end
evaluations.base = Point1_FitnessValue;
evaluations.fhat = fhat_archive;
evaluations.F    = f_archive;
evaluations.count= FEs;
delta.delta1 = delta1;
delta.delta2 = delta2;
end

% function progress(precentage)
% persistent flag = 0;
% str = sprintf('Progress = %d',precentage);
% if flag ~= 1
%     fprintf(1, '%s', str);
%     flag = 1;
% else
%     [char(8)*ones(1,length(str)), str]
%     %fprintf(1, '%s', del);
%     %fprintf(1, '%s', str);
% end
% end

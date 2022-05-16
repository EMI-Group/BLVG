function InitializingPeaks(Dimensions,ii,SeparableWeight)
global MPB;
MPB{ii}.id = ii;
MPB{ii}.Dimension = Dimensions;
if Dimensions>1
    MPB{ii}.PeakNumber = ceil(rand*10);   %%The peaks number between 1 and 10
    MPB{ii}.Weight = 0.5 + rand*1.5;  %% The weight value between 0.5 and 2
else
    MPB{ii}.PeakNumber = ceil(rand*10);
    MPB{ii}.Weight = (0.5 + rand*1.5)*SeparableWeight;   %%Decrease the weight of seperable variables
end
MPB{ii}.ShiftSeverity = 0.5 + rand(MPB{ii}.PeakNumber,1)*2.5;
MPB{ii}.HeightSeverity = 3 + rand(MPB{ii}.PeakNumber,1)*7;
MPB{ii}.WidthSeverity = 0.5 + rand(MPB{ii}.PeakNumber,1)*1;
MPB{ii}.MinCoordinate = -50;
MPB{ii}.MaxCoordinate = 50;
MPB{ii}.MinHeight = 30;
MPB{ii}.MaxHeight = 70;
MPB{ii}.InitialHeight = 50;
MPB{ii}.MinWidth = 1;
MPB{ii}.MaxWidth = 12;
MPB{ii}.InitialWidth = 7;
MPB{ii}.PeaksPosition = MPB{ii}.MinCoordinate + ((MPB{ii}.MaxCoordinate-MPB{ii}.MinCoordinate)*rand(MPB{ii}.PeakNumber,MPB{ii}.Dimension));
MPB{ii}.PeaksHeight = (MPB{ii}.InitialHeight * ones(MPB{ii}.PeakNumber,1)) + (randn(MPB{ii}.PeakNumber,1).*MPB{ii}.HeightSeverity);
MPB{ii}.PeaksWidth = (MPB{ii}.InitialWidth * ones(MPB{ii}.PeakNumber,1)) + (randn(MPB{ii}.PeakNumber,1).*MPB{ii}.WidthSeverity);
end
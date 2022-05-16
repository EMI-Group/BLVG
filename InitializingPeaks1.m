function InitializingPeaks(Dimensions,ii,SeparableWeight)
global MPB;
MPB{ii}.id = ii;
MPB{ii}.Dimension = Dimensions;
MPB{ii}.MinCoordinate      = -50;
MPB{ii}.MaxCoordinate      = 50;
MPB{ii}.MinHeight          = 30;
MPB{ii}.MaxHeight          = 70;
MPB{ii}.MinWidth           = 1;
MPB{ii}.MaxWidth           = 12;
MPB{ii}.PeakNumber         = 10;
MPB{ii}.ShiftSeverity      = ones(MPB{ii}.PeakNumber, 1);
MPB{ii}.HeightSeverity     = 7 * ones(MPB{ii}.PeakNumber, 1);
MPB{ii}.WidthSeverity      = ones(MPB{ii}.PeakNumber, 1);
MPB{ii}.PeaksPosition      = MPB{ii}.MinCoordinate + ((MPB{ii}.MaxCoordinate-MPB{ii}.MinCoordinate) * rand(MPB{ii}.PeakNumber,MPB{ii}.Dimension));
MPB{ii}.PeaksHeight        = MPB{ii}.MinHeight + (MPB{ii}.MaxHeight-MPB{ii}.MinHeight) * rand(MPB{ii}.PeakNumber,1);
MPB{ii}.PeaksWidth         = MPB{ii}.MinWidth + (MPB{ii}.MaxWidth-MPB{ii}.MinWidth) * rand(MPB{ii}.PeakNumber,1);
MPB{ii}.Weight = 1;
end
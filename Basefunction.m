function fit = Basefunction(x,jj)
global  MPB ;
fit = max(MPB{jj}.PeaksHeight -((MPB{jj}.PeaksWidth .* (pdist2(x', MPB{jj}.PeaksPosition))'))) * MPB{jj}.Weight;
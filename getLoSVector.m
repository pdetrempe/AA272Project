function LoS = getLoSVector(x_satellite, x0)
% Use column vectors
LoS = (x_satellite - x0)/norm(x_satellite-x0);
end
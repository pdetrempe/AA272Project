function G = getGeometryMatrix(x_sat_matrix, x0)
% takes in Nx3 satellite matrix and guess at position, returns geometry
% matrix
N = size(x_sat_matrix,1);
G = [];
for i = 1:N
    LoS = getLoSVector(x_sat_matrix(i,:)', x0);
    G = [G; -LoS', 1];
end

end
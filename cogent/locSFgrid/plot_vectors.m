function plot_vectors(block, scale)

dim0 = block(1,1);
dim1 = block(1,2);

n = size(block,1);

R = block(2:n,1);
Z = block(2:n,2);
U = block(2:n,3);
V = block(2:n,4);

dataR = reshape(R, dim0, dim1);
dataZ = reshape(Z, dim0, dim1);
dataU = scale * reshape(U, dim0, dim1);
dataV = scale * reshape(V, dim0, dim1);

quiver(dataR,dataZ,dataU,dataV,0);
axis equal;

end

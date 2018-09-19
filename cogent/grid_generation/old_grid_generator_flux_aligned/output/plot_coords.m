function plot_coords(block)

dim0 = block(1,1)
dim1 = block(1,2)

n = size(block,1);

R = block(2:n,1);
Z = block(2:n,2);

dataR = reshape(R, dim0, dim1);
dataZ = reshape(Z, dim0, dim1);

hold_state = ishold;
hold on

plot_grid(dataR, dataZ, [0 0 0], true, true);
axis equal;

% Restore the hold state
if hold_state == 1
hold on
else
hold off
end

end

function distZ1Star = distance(z1)
global z1Star continuumRad continuumRadVec

continuumRad = 0;
continuumRadVec = zeros(100,1);

% find the distance to the set \A_1
if (norm(z1 - z1Star) <= continuumRad)
    distZ1Star = 0;
else
    distZ1Star = abs(norm(z1) - norm(z1Star + continuumRadVec));
end

end
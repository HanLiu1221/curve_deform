function e = getDeformationEnergy(p1, p2)
%% compute the deformation energe after deforming a curve
e = 0;
for i = 1:length(p1)
    ei = getCurveDeformationEnergy(p1{i}, p2{i});
    e = e + ei;
end
e = e / length(p1);
end

function e = getCurveDeformationEnergy(c1, c2)
%% compute the deformation energe after deforming a curve
e = 0;
for i = 1:length(c1)
    e = e + norm(c1(i, :) - c2(i, :));
end
e = e / length(c1);
end
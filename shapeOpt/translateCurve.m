function deformed = translateCurve(curve, offset)
deformed = zeros(size(curve));
for i = 1 : length(curve)
    deformed(i, :) = curve(i, :) + offset;
end
end
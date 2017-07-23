function drawFeaturePoints(curves, featurIds)
figure;
for i = 1 : length(curves)
    % draw
    plot(curves{i}(:,1), curves{i}(:,2), 'k-');
    hold on;
    plot(curves{i}(featurIds{i}, 1), curves{i}(featurIds{i}, 2), 'r--o', 'LineWidth',2);
    hold on;
    for j = 1:length(featurIds{i})
        text(curves{i}(featurIds{i}(j), 1), curves{i}(featurIds{i}(j), 2), num2str(featurIds{i}(j)));
        hold on;
    end
end
legend('original polyline', 'simplified');
end
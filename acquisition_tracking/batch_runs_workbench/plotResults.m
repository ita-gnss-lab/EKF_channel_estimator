%% Figures configuration
lineWidth = 2;
fontSize = 13;
simulationSteps = 5000; 
epochVector = 1:simulationSteps;

%% Compare different noise levels
delayStd50dB = load('noiseRobustnessTest_50dB\noiseRobustnessTest_50dB_LQGStateRecord.mat');
delayStd38dB = load('noiseRobustnessTest_38dB\noiseRobustnessTest_38dB_LQGStateRecord.mat');
delayStd29dB = load('noiseRobustnessTest_29dB\noiseRobustnessTest_29dB_LQGStateRecord.mat');

figure(Name="STD of the delays", NumberTitle="off");
hold on;
plot(epochVector, ...
    delayStd50dB.LQGStateRecord(1, :), 'LineWidth', lineWidth);
plot(epochVector, ...
    delayStd38dB.LQGStateRecord(1, :), 'LineWidth', lineWidth);
plot(epochVector, ...
    delayStd29dB.LQGStateRecord(1, :), 'LineWidth', lineWidth);
plot(epochVector, 1e-4*ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
legend({"delay STD", "Truth"});
ylabel("Standard deviation of the delays");
xlabel("Epochs (Simulation Steps)");
hold off;

tapsStd50dB = load('noiseRobustnessTest_50dB\noiseRobustnessTest_50dB_channelStateRecord.mat');
tapsStd38dB = load('noiseRobustnessTest_38dB\noiseRobustnessTest_38dB_channelStateRecord.mat');
tapsStd29dB = load('noiseRobustnessTest_29dB\noiseRobustnessTest_29dB_channelStateRecord.mat');
realTaps = load('noiseRobustnessTest_50dB\noiseRobustnessTest_50dB_trueTaps.mat');

channelPlotRows = ceil(sqrt(9));
channelPlotColumns = ceil(9 / channelPlotRows);
figure(Name="Real Channel Tap Estimates 50dB", NumberTitle="off");
tiledlayout(channelPlotRows, channelPlotColumns);
for tapIndex = 10:19
    nexttile;
    hold on;
    plot(epochVector, real(tapsStd50dB.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, real(realTaps.trueTaps(tapIndex - 9)) * ...
        ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
    plot(epochVector, real(tapsStd38dB.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, real(tapsStd29dB.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    hold off;
    title(sprintf("Tap %+d", tapIndex - 9 - 1));
    ylabel("Real");
    xlabel("Epochs");
    set(gca, "FontSize", fontSize);
end
legend({"Estimate", "Truth"});
%% Figures configuration
close all;
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
legend({"50 db-Hz", "38 db-Hz", "29 db-Hz", "Truth"});
ylabel("Delay");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="STD", NumberTitle="off");
hold on;
plot(epochVector, ...
    delayStd50dB.LQGStateRecord(2, :), 'LineWidth', lineWidth);
plot(epochVector, ...
    delayStd38dB.LQGStateRecord(2, :), 'LineWidth', lineWidth);
plot(epochVector, ...
    delayStd29dB.LQGStateRecord(2, :), 'LineWidth', lineWidth);
plot(epochVector, -9.8986558e+05*ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
legend({"50 db-Hz", "38 db-Hz", "29 db-Hz", "Truth"});
ylabel("Phase");
xlabel("Epochs (Simulation Steps)");
hold off;

tapsStd50dB = load('noiseRobustnessTest_50dB\noiseRobustnessTest_50dB_channelStateRecord.mat');
tapsStd38dB = load('noiseRobustnessTest_38dB\noiseRobustnessTest_38dB_channelStateRecord.mat');
tapsStd29dB = load('noiseRobustnessTest_29dB\noiseRobustnessTest_29dB_channelStateRecord.mat');
realTaps = load('noiseRobustnessTest_50dB\noiseRobustnessTest_50dB_trueTaps.mat');

channelPlotRows = ceil(sqrt(4));
channelPlotColumns = ceil(4 / channelPlotRows);
figure(Name="Real Channel Tap Estimates 50dB", NumberTitle="off");
tiledlayout(channelPlotRows, channelPlotColumns);
for tapIndex = 10:13
    nexttile;
    hold on;
    plot(epochVector, real(tapsStd50dB.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, real(tapsStd38dB.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, real(tapsStd29dB.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, real(realTaps.trueTaps(tapIndex - 9)) * ...
        ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
    hold off;
    title(sprintf("Tap -%d", tapIndex - 9 - 1));
    ylabel("Real part of estimated tap");
    xlabel("Epochs");
    set(gca, "FontSize", fontSize);
    legend({"50 db-Hz", "38 db-Hz", "29 db-Hz", "Truth"});
end

channelPlotRows = ceil(sqrt(4));
channelPlotColumns = ceil(4 / channelPlotRows);
figure(Name="Imag Channel Tap Estimates 50dB", NumberTitle="off");
tiledlayout(channelPlotRows, channelPlotColumns);
for tapIndex = 10:13
    nexttile;
    hold on;
    plot(epochVector, imag(tapsStd50dB.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, imag(tapsStd38dB.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, imag(tapsStd29dB.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, imag(realTaps.trueTaps(tapIndex - 9)) * ...
        ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
    hold off;
    title(sprintf("Tap -%d", tapIndex - 9 - 1));
    ylabel("Imaginary part of estimated tap");
    xlabel("Epochs");
    set(gca, "FontSize", fontSize);
    legend({"50 db-Hz", "38 db-Hz", "29 db-Hz", "Truth"});
end

innovationStdRecord50dB = load('noiseRobustnessTest_50dB\noiseRobustnessTest_50dB_innovationStdRecord.mat');
innovationStdRecord38dB = load('noiseRobustnessTest_38dB\noiseRobustnessTest_38dB_innovationStdRecord.mat');
innovationStdRecord29dB = load('noiseRobustnessTest_29dB\noiseRobustnessTest_29dB_innovationStdRecord.mat');
realTaps = load('noiseRobustnessTest_50dB\noiseRobustnessTest_50dB_trueTaps.mat');

figure(Name="STD", NumberTitle="off");
hold on;
plot(epochVector, ...
    innovationStdRecord50dB.innovationStdRecord(:), 'LineWidth', lineWidth);
plot(epochVector, ...
    innovationStdRecord38dB.innovationStdRecord(:), 'LineWidth', lineWidth);
plot(epochVector, ...
    innovationStdRecord29dB.innovationStdRecord(:), 'LineWidth', lineWidth);
plot(epochVector, 0*ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
legend({"50 db-Hz", "38 db-Hz", "29 db-Hz", "Truth"});
ylabel("Innovation STD");
xlabel("Epochs (Simulation Steps)");
hold off;

%% Compare different sampling frequencies 

LQGStateRecord16dB = load('samplesPerChipTest_16dB\samplesPerChipTest_16dB_LQGStateRecord.mat');
LQGStateRecord32dB = load('samplesPerChipTest_32dB\samplesPerChipTest_32dB_LQGStateRecord.mat');
LQGStateRecord64dB = load('samplesPerChipTest_64dB\samplesPerChipTest_64dB_LQGStateRecord.mat');

figure(Name="STys", NumberTitle="off");
hold on;
plot(epochVector, ...
    LQGStateRecord16dB.LQGStateRecord(1, :), 'LineWidth', lineWidth);
plot(epochVector, ...
    LQGStateRecord32dB.LQGStateRecord(1, :), 'LineWidth', lineWidth);
plot(epochVector, ...
    LQGStateRecord64dB.LQGStateRecord(1, :), 'LineWidth', lineWidth);
plot(epochVector, 1e-4*ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
legend({"16 samples", "32 samples", "64 samples", "Truth"});
ylabel("Delay");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="STys", NumberTitle="off");
hold on;
plot(epochVector, ...
    LQGStateRecord16dB.LQGStateRecord(2, :), 'LineWidth', lineWidth);
plot(epochVector, ...
    LQGStateRecord32dB.LQGStateRecord(2, :), 'LineWidth', lineWidth);
plot(epochVector, ...
    LQGStateRecord64dB.LQGStateRecord(2, :), 'LineWidth', lineWidth);
plot(epochVector, -9.8986558e+05 *ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
legend({"16 samples", "32 samples", "64 samples", "Truth"});
ylabel("Delay");
xlabel("Epochs (Simulation Steps)");
hold off;

tapsStd16 = load('samplesPerChipTest_16dB\samplesPerChipTest_16dB_channelStateRecord.mat');
tapsStd32 = load('samplesPerChipTest_32dB\samplesPerChipTest_32dB_channelStateRecord.mat');
tapsStd64 = load('samplesPerChipTest_64dB\samplesPerChipTest_64dB_channelStateRecord.mat');
realTaps = load('samplesPerChipTest_16dB\samplesPerChipTest_16dB_trueTaps.mat');

channelPlotRows = ceil(sqrt(4));
channelPlotColumns = ceil(4 / channelPlotRows);
figure(Name="Real Channel Tap different sampling", NumberTitle="off");
tiledlayout(channelPlotRows, channelPlotColumns);
for tapIndex = 10:13
    nexttile;
    hold on;
    plot(epochVector, real(tapsStd16.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, real(tapsStd32.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, real(tapsStd64.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, real(realTaps.trueTaps(tapIndex - 9)) * ...
        ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
    hold off;
    title(sprintf("Tap -%d", tapIndex - 9 - 1));
    ylabel("Real part of estimated tap");
    xlabel("Epochs");
    set(gca, "FontSize", fontSize);
    legend({"16 samples", "32 samples", "64 samples", "Truth"});
end

channelPlotRows = ceil(sqrt(4));
channelPlotColumns = ceil(4 / channelPlotRows);
figure(Name="Imag Channel Tap different sampling", NumberTitle="off");
tiledlayout(channelPlotRows, channelPlotColumns);
for tapIndex = 10:13
    nexttile;
    hold on;
    plot(epochVector, imag(tapsStd16.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, imag(tapsStd32.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, imag(tapsStd64.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, imag(realTaps.trueTaps(tapIndex - 9)) * ...
        ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
    hold off;
    title(sprintf("Tap -%d", tapIndex - 9 - 1));
    ylabel("Imaginary part of estimated tap");
    xlabel("Epochs");
    set(gca, "FontSize", fontSize);
    legend({"50 db-Hz", "38 db-Hz", "29 db-Hz", "Truth"});
end


channelPlotRows = ceil(sqrt(4));
channelPlotColumns = ceil(4 / channelPlotRows);
figure(Name="Real Channel Tap different sampling", NumberTitle="off");
tiledlayout(channelPlotRows, channelPlotColumns);
for tapIndex = 10:13
    nexttile;
    hold on;

    plot(real(tapsStd16.channelStateRecord(tapIndex, end)), ...
    imag(tapsStd16.channelStateRecord(tapIndex, end)), ...
    'o', 'LineWidth', lineWidth, 'MarkerSize', 8);

    text(real(tapsStd16.channelStateRecord(tapIndex, end)), ...
    imag(tapsStd16.channelStateRecord(tapIndex, end)), ...
    sprintf(" -%d", tapIndex - q - 1));

    plot(real(tapsStd32.channelStateRecord(tapIndex, end)), ...
    imag(tapsStd32.channelStateRecord(tapIndex, end)), ...
    'o', 'LineWidth', lineWidth, 'MarkerSize', 8);

    text(real(tapsStd32.channelStateRecord(tapIndex, end)), ...
    imag(tapsStd32.channelStateRecord(tapIndex, end)), ...
    sprintf(" -%d", tapIndex - q - 1));

    plot(real(tapsStd64.channelStateRecord(tapIndex, end)), ...
    imag(tapsStd64.channelStateRecord(tapIndex, end)), ...
    'o', 'LineWidth', lineWidth, 'MarkerSize', 8);

    text(real(tapsStd64.channelStateRecord(tapIndex, end)), ...
    imag(tapsStd64.channelStateRecord(tapIndex, end)), ...
    sprintf(" -%d", tapIndex - q - 1));

    plot(real(realTaps.trueTaps(tapIndex - 9)), ...
    imag(realTaps.trueTaps(tapIndex - 9)), 'x', ...
    'LineWidth', lineWidth, 'MarkerSize', 10);
     
    text(real(realTaps.trueTaps(tapIndex - 9)), ...
    imag(realTaps.trueTaps(tapIndex - 9)), ...
    sprintf(" -%d", tapIndex - q - 1));
    buffer = 0.05 * 0.8^(tapIndex - 9);
    dX = max(abs(real(realTaps.trueTaps(tapIndex - 9)) - xlim));
    dY = max(abs(imag(realTaps.trueTaps(tapIndex - 9)) - ylim));
    xlim(real(realTaps.trueTaps(tapIndex - 9))+[(-dX-buffer) (dX+buffer)]);
    ylim(imag(realTaps.trueTaps(tapIndex - 9))+[(-dY-buffer) (dY+buffer)]);

    hold off;
    title(sprintf("Tap -%d", tapIndex - 9 - 1));
    grid on;
    xlabel("Real");
    ylabel("Imaginary");
    set(gca, "FontSize", fontSize);
    legend({"16 samples", "32 samples", "64 samples", "Truth"});
end

figure(Name="Real Channel Tap different sampling", NumberTitle="off");
hold on;
for tapIndex = 10:18
    plot(real(realTaps.trueTaps(tapIndex - 9)), ...
    imag(realTaps.trueTaps(tapIndex - 9)), 'x', ...
    'LineWidth', lineWidth, 'MarkerSize', 10, Color=[0.5, 0.2, 0.6]);
    grid on;
    xlabel("Real");
    ylabel("Imaginary");
    set(gca, "FontSize", fontSize);
end
legend({"Simulation Taps"});
hold off;

%% Compare different number of correltaors

LQGStateRecord9dB = load('numberOfCorrelatorsTest_9dB\numberOfCorrelatorsTest_9dB_LQGStateRecord.mat');
LQGStateRecord15dB = load('numberOfCorrelatorsTest_15dB\numberOfCorrelatorsTest_15dB_LQGStateRecord.mat');
LQGStateRecord21dB = load('numberOfCorrelatorsTest_21dB\numberOfCorrelatorsTest_21dB_LQGStateRecord.mat');

figure(Name="STys", NumberTitle="off");
hold on;
plot(epochVector, ...
    LQGStateRecord9dB.LQGStateRecord(1, :), 'LineWidth', lineWidth);
plot(epochVector, ...
    LQGStateRecord15dB.LQGStateRecord(1, :), 'LineWidth', lineWidth);
plot(epochVector, ...
    LQGStateRecord21dB.LQGStateRecord(1, :), 'LineWidth', lineWidth);
plot(epochVector, 1e-4*ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
legend({"9 correlators", "15 correlators", "21 correlators", "Truth"});
ylabel("Delay");
xlabel("Epochs (Simulation Steps)");
hold off;

figure(Name="STys", NumberTitle="off");
hold on;
plot(epochVector, ...
    LQGStateRecord9dB.LQGStateRecord(2, :), 'LineWidth', lineWidth);
plot(epochVector, ...
    LQGStateRecord15dB.LQGStateRecord(2, :), 'LineWidth', lineWidth);
plot(epochVector, ...
    LQGStateRecord21dB.LQGStateRecord(2, :), 'LineWidth', lineWidth);
plot(epochVector, -9.8986558e+05 *ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
legend({"9 correlators", "15 correlators", "21 correlators", "Truth"});
ylabel("Delay");
xlabel("Epochs (Simulation Steps)");
hold off;

tapsStd9dB = load('numberOfCorrelatorsTest_9dB\numberOfCorrelatorsTest_9dB_channelStateRecord.mat');
tapsStd15dB = load('numberOfCorrelatorsTest_15dB\numberOfCorrelatorsTest_15dB_channelStateRecord.mat');
tapsStd21dB = load('numberOfCorrelatorsTest_21dB\numberOfCorrelatorsTest_21dB_channelStateRecord.mat');
realTaps = load('samplesPerChipTest_16dB\samplesPerChipTest_16dB_trueTaps.mat');

channelPlotRows = ceil(sqrt(4));
channelPlotColumns = ceil(4 / channelPlotRows);
figure(Name="g", NumberTitle="off");
tiledlayout(channelPlotRows, channelPlotColumns);
for tapIndex = 10:13
    nexttile;
    hold on;
    plot(epochVector, real(tapsStd9dB.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, real(tapsStd15dB.channelStateRecord(tapIndex + 6, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, real(tapsStd21dB.channelStateRecord(tapIndex + 12, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, real(realTaps.trueTaps(tapIndex - 9)) * ...
        ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
    hold off;
    title(sprintf("Tap -%d", tapIndex - 9 - 1));
    ylabel("Real part of estimated tap");
    xlabel("Epochs");
    set(gca, "FontSize", fontSize);
    legend({"9 correlators", "15 correlators", "21 correlators", "Truth"});
end

channelPlotRows = ceil(sqrt(4));
channelPlotColumns = ceil(4 / channelPlotRows);
figure(Name="ng", NumberTitle="off");
tiledlayout(channelPlotRows, channelPlotColumns);
for tapIndex = 10:13
    nexttile;
    hold on;
    plot(epochVector, imag(tapsStd9dB.channelStateRecord(tapIndex, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, imag(tapsStd15dB.channelStateRecord(tapIndex +6, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, imag(tapsStd21dB.channelStateRecord(tapIndex +12, :)), ...
        'LineWidth', lineWidth);
    plot(epochVector, imag(realTaps.trueTaps(tapIndex - 9)) * ...
        ones(1, simulationSteps), '--', 'LineWidth', lineWidth);
    hold off;
    title(sprintf("Tap -%d", tapIndex - 9 - 1));
    ylabel("Imaginary part of estimated tap");
    xlabel("Epochs");
    set(gca, "FontSize", fontSize);
    legend({"9 correlators", "15 correlators", "21 correlators", "Truth"});
end

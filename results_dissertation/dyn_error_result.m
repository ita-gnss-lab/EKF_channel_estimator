clearvars; clc; close all;

addpath(genpath(fullfile("..", "..","EKF_channel_estimator")));

load default_config.mat
exactEstimate = [1e-4, ...
    configuration.dopplerProfile(1), ...
    2*pi*(configuration.dopplerProfile(2)), ...
    2*pi*configuration.dopplerProfile(3)].';

configuration.steps = 1000;
delayRecord = zeros(10,configuration.steps);
delayTrue = zeros(10,configuration.steps);
configuration.dopplerProfile(3) = 500;
for i = 1:5
    rng(2643726);
    initialEstimate = exactEstimate;
    initialEstimate(4) = 2*pi*500*(1 + 0.1*i);
    [delayRecord(i,:), delayTrue(i,:)] = lqg_function(configuration, initialEstimate);
end
%%
window = 30;
errorRMSE= zeros(10, configuration.steps - window);
for i = 1:5
    for k = 1:(configuration.steps - window)
        errorRMSE(i, k) = rmse(delayRecord(i,k:(k+window)), delayTrue(i,k:(k+window))) * (100/1e-4);
    end
end

lineWidth = 2;
fontSize = 13;
figure(Name="Delay RMSE", NumberTitle="off");
hold on;
for i = 1:5
    plot(errorRMSE(i, :), 'LineWidth', lineWidth);
end
legend({"10%", "20%", "30%", "40%", "50%"}, 'Location', 'northwest');
ylabel("Delay RMSE in %");
ylim([0 0.25]);
xlabel("Epochs (Simulation Steps)");
hold off;

function [LOS_phase, LOS_delay, time] = get_LOS_dynamics(time, doppler_profile, carrier_freq)
%GET_LOS_DYNAMICS Summary of this function goes here
%   Detailed explanation goes here

initial_condition = doppler_profile(1) * ones(size(time))';
% LOS phase's Taylor's expansion given the doppler derivatives
if length(doppler_profile) > 1
    orders = 1 : (length(doppler_profile) - 1);
    taylors_factors = (time'.^orders)./factorial(orders);
    taylors_expansion = 2*pi*sum(doppler_profile(orders + 1) .* taylors_factors, 2);
    LOS_phase = initial_condition + taylors_expansion;  
else
    LOS_phase = initial_condition;
end

% LOS delay from carrier phase
LOS_delay = -LOS_phase ./ (2*pi*carrier_freq);

end


function [derivative] = autocorrelationDerivative(delay, Chipping_Frequency, Sampling_Frequency)
%AUTOCORRELATIONDERIVATIVE Summary of this function goes here
%   Detailed explanation goes here
    Band_Width = Sampling_Frequency/2;
    Integral_Interval = Band_Width * (-1: 1/1000 : 1);
    Integral_Interval(1001) = (1e-18); 
    
    SIN2 = sin(pi/Chipping_Frequency * Integral_Interval).^2;
    MULT = SIN2./(pi*Integral_Interval/Chipping_Frequency);
    SIN  = sin(2*pi * delay * Integral_Interval);
    Function = 2*MULT.*SIN;
    derivative = -trapz(Integral_Interval, Function);
end


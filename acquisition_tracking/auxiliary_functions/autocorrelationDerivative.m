function [derivative] = autocorrelationDerivative(inputArg1,inputArg2)
%AUTOCORRELATIONDERIVATIVE Summary of this function goes here
%   Detailed explanation goes here
    Band_Width = Sampling_Frequency/2;
    Integral_Interval = -Band_Width: Band_Width/1000 : Band_Width;
    Integral_Interval(1001) = (1e-18);
    
    SIN2 = sin(pi/Chipping_Frequency * Integral_Interval).^2;
    MULT = SIN2./(pi*Integral_Interval/Chipping_Frequency);
    SIN  = sin(2*pi/Chipping_Frequency * Integral_Interval * (j-i)*Maximum_Spacement/Number_of_Spacings);
    Function = 2*MULT.*SIN;
    derivative = -trapz(Integral_Interval, Function);
end


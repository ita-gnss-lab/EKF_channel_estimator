function [signal, nextSampleIndex, File_Ended] = get_from_DAT(fileName, inputSampleIndex, samplesTotal)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[signal, outputSampleIndex, File_Ended] = read_gr_complex_binary ( ...
    fileName, ...
    inputSampleIndex, ...
    samplesTotal);

if File_Ended == true
    return;
end
disp(inputSampleIndex);

nextSampleIndex = outputSampleIndex;

end


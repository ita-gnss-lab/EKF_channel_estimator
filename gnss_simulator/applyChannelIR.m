function y = applyChannelIR(x, h)
% APPLYCHANNELIR  Offline FFT convolution through an integer-sample FIR channel.
%   y = applyChannelIR(x, h)
%
% Description
%   Filters x using the impulse response h with fftfilt() (overlap-add).
%   Intended for recorded GNSS baseband where tap delays are integer
%   multiples of Ts and encoded directly in h.
%
% Inputs
%   x : vector (real or complex) — input baseband signal
%   h : vector (real or complex) — FIR where h(k) is the tap at delay k-1
%
% Output
%   y : column vector — filtered signal (length ≈ length(x)+length(h)-1)
%
% Fixed behavior
%   - Uses the provided impulse response without tapering or renormalization.
%   - Uses fftfilt() for efficient O(N log N) convolution.

    x = x(:);
    h = h(:);

    % FFT-based convolution only
    y = fftfilt(h, x);
end

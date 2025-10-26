function y = applyChannelIR(x, h)
% APPLYCHANNELIR  Offline FFT convolution through an integer-sample FIR channel.
%   y = applyChannelIR(x, h)
%
% Description
%   Applies a Hamming window to the impulse response h, renormalizes it to
%   preserve the L2 (energy) of the original taps, and filters x using
%   fftfilt() (overlap-add). Intended for recorded GNSS baseband where tap
%   delays are integer multiples of Ts and encoded directly in h.
%
% Inputs
%   x : vector (real or complex) — input baseband signal
%   h : vector (real or complex) — FIR where h(k) is the tap at delay k-1
%
% Output
%   y : column vector — filtered signal (length ≈ length(x)+length(h)-1)
%
% Fixed behavior
%   - Uses hamming(length(h)) to taper the IR.
%   - Renormalizes to keep sum(|h|^2) unchanged (L2 energy).
%   - Uses fftfilt() for efficient O(N log N) convolution.

    x = x(:);
    h = h(:);

    % Hamming window
    w  = hamming(numel(h));
    hW = h .* w;

    % L2 (energy) preservation
    e0 = sum(abs(h).^2);
    eW = sum(abs(hW).^2);
    hW = hW * sqrt(e0 / eW);

    % FFT-based convolution only
    y = fftfilt(hW, x);
end

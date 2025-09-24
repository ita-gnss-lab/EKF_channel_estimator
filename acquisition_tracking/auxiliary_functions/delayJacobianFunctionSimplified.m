function [hTau, dphiLm, deltaLm] = delayJacobianFunctionSimplified( ...
        epsTauHat, hHatL, Ts, q, Tc, scale)
% delayJacobianFunctionSimplified
% Computes the delay Jacobian vector hTau (size (2q+1)x1) whose m-th entry is
%   hTau(m) = scale * sum_{l=0}^L hHatL(l+1) * dPhi_pp(deltaLm(l+1,m))
% with
%   deltaLm(l+1,m) = epsTauHat + (l + m) * Ts,
%   dPhi_pp(Δ) = -(1/Tc)*sign(Δ),  for 0 < |Δ| < Tc; and 0 elsewhere.
%
% Inputs
%   epsTauHat : scalar delay error estimate ε̂_τ[k|k−1]
%   hHatL     : (L+1)x1 vector of channel coeffs for l = 0..L (complex ok)
%   Ts        : sampling period
%   q         : half-window for m → m ∈ {-q, …, q}
%   Tc        : (optional) chip duration; default 1 (Δ already normalized)
%   scale     : (optional) global complex scale; default 1
%               e.g., scale = sqrt(2*Ps) * exp(1j*epsPhiHat)
%
% Outputs
%   hTau    : (2q+1)x1 Jacobian vector over m
%   dphiLm  : (L+1)x(2q+1) matrix with dPhi_pp(deltaLm)
%   deltaLm : (L+1)x(2q+1) matrix with Δ_{l,m}
%
% Notes
% - At Δ ∈ {-Tc, 0, Tc} the derivative is set to 0 (subgradient choice).
% - Fully vectorized; requires implicit expansion (R2016b+).

    if nargin < 5 || isempty(Tc),    Tc = 1;   end
    if nargin < 6 || isempty(scale), scale = 1; end

    % Indices
    L = numel(hHatL) - 1;
    l = (0:L).';       % (L+1)x1
    m = (-q:q);        % 1x(2q+1)

    % Δ_{l,m} = ε̂_τ + (l + m)Ts
    deltaLm = epsTauHat + (l + m) * Ts;

    % dΦ_pp(Δ): derivative of triangular autocorrelation
    dphiLm = zeros(size(deltaLm), 'like', deltaLm);
    interior = (abs(deltaLm) < Tc) & (deltaLm ~= 0);
    dphiLm(interior) = -(1./Tc) .* sign(deltaLm(interior));
    % elsewhere remains 0

    % Per-(l,m) contributions and sum over l for each m
    contribLm = dphiLm .* hHatL(:);   % (L+1)x(2q+1)
    hTau = scale .* sum(contribLm, 1).';  % (2q+1)x1
end

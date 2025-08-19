% get_wiener_model_symbolic
%
% Syntax:
%   get_wiener_model_symbolic()
%
% Description:
%   Computes the symbolic version of the discrete Wiener state space model.
%   All equations used here are referenced to [1].
%   
% Inputs:
%   - No Inputs
% Outputs:
%   - It outputs the state transition matrix F_W, and the discrete noise 
%   covariance matrices that composes Q_D.
% 
% Notes:
%   - In [1], the frequency shift and drift are modeled in Hz, while here
%   we modeled them in rad/s and rad/s², respectively. That is why the 2\pi
%   factor does not exist in this code.
%
% References:
%   [1] R. V. Pacelli, R. D. L. Florindo, F. Antreich, and A. M. P. 
%       De Lucena, “An all-digital coherent AFSK demodulator for CubeSat 
%       applications,” Digital Signal Processing, vol. 162, p. 105147, 
%       Jul. 2025, doi: 10.1016/j.dsp.2025.105147.
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

% Defining symbolic variables
syms beta_delay T tau real

% SEE: Fw in [1, Eq. 30].
% NOTE: Now we also model the delay, so that is why it is a littlbe bit
% different.
F_w_continuous_time = sym([0 0 beta_delay 0;
         0 0 1    0;
         0 0 0    1;
         0 0 0    0]);

% Discrete-time transition matrix
% SEE: [1, Eq. 33], for t - tau = T.
F_W = expm(F_w_continuous_time * T);

% Display F_W in the command line
disp(F_W);

% SEE: [1, Eq. 33].
A_T_minus_tau = expm(F_w_continuous_time * (T - tau));

% Note that we can break Q_xi as a sum of 4 matrix in this case. I'm
% not putting the variances {\sigma_1, \sigma_2, ...}, in order to
% avoid visual cluster in the plots.
Q_xi_1 = diag([1,0,0,0]);
Q_xi_2 = diag([0,1,0,0]);
Q_xi_3 = diag([0,0,1,0]);
Q_xi_4 = diag([0,0,0,1]);

% SEE: [1, Eq. 38].
QD_1 = int(A_T_minus_tau * Q_xi_1 * A_T_minus_tau.', tau, 0, T);
QD_2 = int(A_T_minus_tau * Q_xi_2 * A_T_minus_tau.', tau, 0, T);
QD_3 = int(A_T_minus_tau * Q_xi_3 * A_T_minus_tau.', tau, 0, T);
QD_4 = int(A_T_minus_tau * Q_xi_4 * A_T_minus_tau.', tau, 0, T);

disp(QD_1);
disp(QD_2);
disp(QD_3);
disp(QD_4);
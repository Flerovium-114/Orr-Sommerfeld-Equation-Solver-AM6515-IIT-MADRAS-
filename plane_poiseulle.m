% Plane Poiseuille flow Orr-Sommerfeld solver
clear; clc; close all;

% Parameters
Re = 10000;          % Reynolds number
n = 100;             % Number of Chebyshev collocation points
alpha_range = linspace(-0.1, 0.1, 100); % Range of alpha (wavenumbers)
Re_range = linspace(2000, 10000, 20);  % Range of Reynolds numbers

% Chebyshev differentiation matrix
[D1, x] = chebyshev_diff_matrix(n);
D2 = D1^2;           % Second derivative
D2 = D2(2:end-1, 2:end-1); % Remove boundary rows for v'' BCs
I = eye(n-1);        % Identity matrix for inner points

% Base velocity profile for plane Poiseuille flow
U = diag(1 - x.^2);      % U(y) = 1 - y^2
U = U(2:end-1, 2:end-1); % Remove boundary rows/columns
U_dd = diag(-2 * ones(n-1, 1)); % U'' = -2 (constant for Poiseuille flow)

% Preallocate for storing results
growth_rates = zeros(1, length(alpha_range));
phase_speeds = zeros(1, length(alpha_range));
eigReal = [];
eigImag = [];
neutral_Re = [];

% Loop over wavenumbers
for j = 1:length(alpha_range)
    alpha = alpha_range(j);

    % Orr-Sommerfeld matrices
    L = -1i * alpha * Re * (D2 - alpha^2);  % Orr-Sommerfeld operator
    M = -1i * alpha * Re * (U * (D2 - alpha^2) - U_dd) + (D2 - alpha^2)^2; % RHS operator

    % Apply boundary conditions for v and v'
    A = [L, zeros(size(L)); 
         zeros(size(L)), I];
    B = [M, zeros(size(M)); 
         zeros(size(M)), I];

    % Solve generalized eigenvalue problem
    [eigVec, eigVal] = eig(A, B);
    eigVals = diag(eigVal);

    % Extract and store real and imaginary parts of eigenvalues
    eigReal = [eigReal; real(eigVals)];
    eigImag = [eigImag; imag(eigVals)];
    

    % Extract growth rate (imaginary part of eigenvalues)
    [maxGrowth, maxIdx] = max(imag(eigVals));
    growth_rates(j) = maxGrowth;
    phase_speeds(j) = real(eigVals(maxIdx)) / alpha;
end

% Plot results
figure;
plot(alpha_range, growth_rates, '.', 'LineWidth', 1.5);
xlabel('\alpha'); ylabel('Growth rate (Imag(\omega))');
title('Stability Curve for Plane Poiseuille Flow');
grid on;

figure;
plot(alpha_range, phase_speeds, '.', 'LineWidth', 1.5);
xlabel('\alpha'); ylabel('Phase speed (Real(\omega))');
grid on;

figure;
plot(phase_speeds, growth_rates, '.', 'LineWidth', 1.5);
xlabel('Phase Speed (Real part of \omega))'); ylabel('Growth Rate (Imag(\omega))');
grid on;

% Supporting function: Chebyshev differentiation matrix
function [D, x] = chebyshev_diff_matrix(n)
    % Returns the Chebyshev differentiation matrix D and grid points x
    x = cos(pi * (0:n) / n)'; % Chebyshev points
    c = [2; ones(n-1, 1); 2] .* (-1).^(0:n)'; % Scaling factors
    X = repmat(x, 1, n+1); % Repeated x-grid
    dX = X - X'; % Difference matrix
    D = (c * (1 ./ c)') ./ (dX + eye(n+1)); % Differentiation matrix
    D = D - diag(sum(D')); % Set diagonal elements
end

function [S, ST] = gif1_sensitivity_analysis(M, odes)

%% Parameters
% QC
k_1_qc = 0.01 + (3000 - 0.01) * rand(2 * M, 1);
k_2_qc = -7+(7+7)*rand(2*M,1);
k_3_qc = 0.01 + (3000 - 0.01) * rand(2 * M, 1);
k_8_qc = 0.05 + (1000 - 0.05) * rand(2 * M, 1);
K_D2_qc = 1.5 + (1000 - 1.5) * rand(2 * M, 1);
K_D3_qc = 1.5 + (1000 - 1.5) * rand(2 * M, 1);
K_D4_qc = 1.5 + (1000 - 1.5) * rand(2 * M, 1);
d_3_qc = 0.01 + (1000 - 0.01) * rand(2 * M, 1);
d_8_qc = 0.01 + (200 - 0.01) * rand(2 * M, 1);
d_12_qc = 0.01 + (200 - 0.01) * rand(2 * M, 1);

% CEI
k_2_cei = -7+(7+7)*rand(2*M,1);
k_3_cei = 0.01 + (3000 - 0.01) * rand(2 * M, 1);
k_5_cei = 0.01 + (3000 - 0.01) * rand(2 * M, 1);
k_8_cei = 0.05 + (1000 - 0.05) * rand(2 * M, 1);
K_D2_cei = 1.5 + (1000 - 1.5) * rand(2 * M, 1);
K_D3_cei = 1.5 + (1000 - 1.5) * rand(2 * M, 1);
K_D4_cei = 1.5 + (1000 - 1.5) * rand(2 * M, 1);
d_3_cei = 0.01 + (1000 - 0.01) * rand(2 * M, 1);
d_5_cei = 0.01 + (1000 - 0.01) * rand(2 * M, 1);
d_8_cei = 0.01 + (200 - 0.01) * rand(2 * M, 1);
d_11_cei = 0.01 + (1000 - 0.01) * rand(2 * M, 1);
d_12_cei = 0.01 + (1000 - 0.01) * rand(2 * M, 1);

% VASC
k_4_vasc = 0.05 + (1000 - 0.05) * rand(2 * M, 1);
K_D1_vasc = 1.5 + (1000 - 1.5) * rand(2 * M, 1);
d_4_vasc = 0.01 + (1000 - 0.01) * rand(2 * M, 1);

% ENDO
k_2_endo = -7+(7+7)*rand(2*M,1);
k_3_endo = 0.01 + (3000 - 0.01) * rand(2 * M, 1);
k_8_endo = 0.01 + (3000 - 0.01) * rand(2 * M, 1);
K_D2_endo = 1.5 + (1000 - 1.5) * rand(2 * M, 1);
K_D3_endo = 1.5 + (1000 - 1.5) * rand(2 * M, 1);
K_D4_endo = 1.5 + (1000 - 1.5) * rand(2 * M, 1);
d_3_endo = 0.01 + (1000 - 0.01) * rand(2 * M, 1);
d_8_endo = 0.01 + (1000 - 0.01) * rand(2 * M, 1);
d_12_endo = 0.01 + (200 - 0.01) * rand(2 * M, 1);

params = [k_1_qc,k_2_qc,k_3_qc,k_8_qc,K_D2_qc,K_D3_qc,K_D4_qc,d_3_qc,d_8_qc,d_12_qc...
    k_2_cei,k_3_cei,k_5_cei,k_8_cei,K_D2_cei,K_D3_cei,K_D4_cei,d_3_cei,d_5_cei...
    d_8_cei,d_11_cei,d_12_cei,...
    k_4_vasc,K_D1_vasc,d_4_vasc,...
    k_2_endo,k_3_endo,k_8_endo,K_D2_endo,K_D3_endo,K_D4_endo,d_3_endo,d_8_endo,...
    d_12_endo];
[~, num_params] = size(params);
A = params(1:M, :);
B = params((M + 1):(2 * M), :);

% Preallocate
ya = zeros(M, odes);
yb = zeros(M, odes);
yc = zeros(M, odes);
S = zeros(size(A, 2), odes);
ST = zeros(size(A, 2), odes);

%% Calculate Sobol indices for each parameter
% Comp time
tic;

for i = 1:size(A, 2)
    % Create C, which is equal to A except that the ith column of B is
    % replaced by the ith column of B
    C = A;
    C(:, i) = B(:, i);
    % For each model simulation
    for j = 1:M
        % Save the parameters for this simulation and feed them into the
        % ODE model
        aparams = A(j, :);
        bparams = B(j, :);
        cparams = C(j, :);
        myfuna = @(t, x) gif1_dy(t, x, aparams(1:num_params));
        myfunb = @(t, x) gif1_dy(t, x, bparams(1:num_params));
        myfunc = @(t, x) gif1_dy(t, x, cparams(1:num_params));
        % Use ODE45 to solve the ODE
        y0a = [161.41,107.18,307.49,127.62,0,...
            161.41,107.18,307.49,127.62,0,7.16,...
            161.41,107.18,...
            107.18,307.49,127.62,0];
        y0b = [161.41,107.18,307.49,127.62,0,...
            161.41,107.18,307.49,127.62,0,7.16,...
            161.41,107.18,...
            107.18,307.49,127.62,0];
        y0c = [161.41,107.18,307.49,127.62,0,...
            161.41,107.18,307.49,127.62,0,7.16,...
            161.41,107.18,...
            107.18,307.49,127.62,0];
        [T, YA] = ode45(myfuna, 0:(1 / 60):(24 - (1 / 60)), y0a, []);
        [~, YB] = ode45(myfunb, 0:(1 / 60):(24 - (1 / 60)), y0b, []);
        [~, YC] = ode45(myfunc, 0:(1 / 60):(24 - (1 / 60)), y0c, []);
        % Use trapz to numerically integrate the solution over time. This
        % is used to calcualte the Sobol index. trapz integrates each
        % column separately
        % Sometimes, the solver might not make it all the way through the
        % time series (this is rare). In this case, we will just use the
        % values from a random previous iteration
        if size(YA, 1) ~= length(T) || size(YB, 1) ~= length(T)...
                || size(YC, 1) ~= length(T)
            iter = floor(1 + rand*(j-2));
            ya(j, :) = ya(iter, :);
            yb(j, :) = yb(iter, :);
            yc(j, :) = yc(iter, :);
        else
            ya(j, :) = trapz(T, YA);
            yb(j, :) = trapz(T, YB);
            yc(j, :) = trapz(T, YC);
        end
    end
    % Calculate sobol indices for each of the responses
    % Estimators from Saltelli et al, 2010
    for k = 1:odes
        amc = (ya(:, k) - yc(:, k))' * (ya(:, k) - yc(:, k));
        S(i, k) = (1 / M * yb(:, k)' * (yc(:, k) - ya(:, k))) / var(yb(:, k));
        ST(i, k) = amc / (2 * M) / var(yb(:, k));
    end
end

% Comp time
toc;

end


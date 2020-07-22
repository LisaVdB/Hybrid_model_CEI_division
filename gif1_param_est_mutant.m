% Function to run parameter estimtion

function [q_array, fvals] = gif1_param_est_mutant(M)
%% Initialization
% Use latin hypercube for parameter estimation
% LHS provides faster convergence than Monte Carlo
mymat = lhsdesign(M, 11);

% Build matrix of initial parameter estimates
% Ranges based on initial estimates that approximate the data
k_3_endo = 1000 + (17000 - 1000) * mymat(:, 1);
d_3_endo = 0.1 + (15 - 0.1) * mymat(:, 2);
d_5_cei = 1 + (15 - 1) * mymat(:, 3);
k_6_cei = 100 + (6000 - 100) * mymat(:, 4);
d_6_cei = 0.1 + (150 - 0.1) * mymat(:, 5);
k_3_cei = 1 + (600 - 1) * mymat(:, 1);

q_initial = [k_3_endo, d_3_endo, d_5_cei, k_6_cei, d_6_cei, k_3_cei];

% Upper and lower bounds of parameters 
ub = max(q_initial);
lb = min(q_initial);

% Build arrays to store results
q_array = zeros(M, size(q_initial, 2));
fvals = zeros(M, 1);

%% Fit square residual using simulated annealing

% Start comp time
tic;

for i = 1:M
    fprintf("\nRunning simulation %d of %d...\n", i, M);
    % Stop after max_time seconds or max_calc calculations
    max_time = 360;
    max_calc = 39000;
    options = saoptimset('TimeLimit', max_time); %#ok<*SAOPT>
    options.MaxFunctionEvaluations = max_calc;
    [q, fval, ~, ~] = simulannealbnd(@residual, q_initial(i, :), lb, ub, options);
    % Save results
    q_array(i, :) = q;
    fvals(i) = fval;
end

fprintf("Done. ");
% End comp time
toc;

%% Residual calculation function
    function R = residual(q)
        % Evaluate ODE
        % 4D to 4D 8H
        y0 = [161.41, 0, 135.41, 127.62, 0,...
            0, 0, 307.49, 45.93, 0, 7.16,...
            0, 107.18,...
            0, 87.80, 127.62, 0];
        q0 = [q, 14.5474, -.00575195]; % Adjust params at each time break
        myfun = @(t, x)gif1_dy_est_mutant(t, x, q0);
        [T0, Y0] = ode45(myfun, 0:(1/60):(8-(1/60)), y0, []);
        
        % 4D 8H to 4D 16H
        y01 = Y0(size(Y0, 1), :);
        q1 = [q, 14.5808, -0.0615631];
        myfun = @(t, x)gif1_dy_est_mutant(t, x, q1);
        [T1, Y1] = ode45(myfun, 8:(1/60):(16-(1/60)), y01, []);
        
        % 4D 16H to 5D
        y02 = Y1(size(Y1, 1), :);
        q2 = [q, 14.8001, 0.010347];
        myfun = @(t, x)gif1_dy_est_mutant(t, x, q2);
        [T2, Y2] = ode45(myfun, 16:(1/60):(24-(1/60)), y02, []);
        
        % 5D - 5D 8H
        y03 = Y2(size(Y2, 1), :);
        q3 = [q, 14.5075, 0.0548233];
        myfun = @(t, x)gif1_dy_est_mutant(t, x, q3);
        [T3, Y3] = ode45(myfun, 24:(1/60):(32-(1/60)), y03, []);
        
        % 5D 8H - 5D 16H
        y04 = Y3(size(Y3, 1), :);
        q4 = [q, 14.7837, 0.0775361];
        myfun = @(t,x)gif1_dy_est_mutant(t, x, q4);
        [T4, Y4] = ode45(myfun, 32:(1/60):(40-(1/60)), y04, []);
        
        % 5D 16H - 6D
        y05 = Y4(size(Y4, 1), :);
        q5 = [q, 14.5807, -0.1128];
        myfun = @(t,x)gif1_dy_est_mutant(t, x, q5);
        [T5, Y5] = ode45(myfun, 40:(1/60):(48-(1/60)), y05, []);
        
        % Combine results
        Tf = [T0; T1; T2; T3; T4; T5];
        Yf = [Y0; Y1; Y2; Y3; Y4; Y5];
        % Get solutions
        WOX5 = Yf(:, 1);
        SHR = Yf(:, 13);
        AN3_cei = Yf(:, 8);
        AN3_qc = Yf(:, 3);
        AN3_endo = Yf(:, 15);
        SCR_endo = Yf(:, 16);
        SCR_cei = Yf(:, 9);
        CYCD6 = Yf(:, 11);
        X = Yf(:, 6);
        
        % Calculate residuals (2880 = 48 H - one timestep)
        % R1: WOX5 time course data
        R1 = (WOX5(Tf==8)-80.47)^2 + (WOX5(Tf==16)-52.38)^2 +...
            (WOX5(Tf==24)-196.32)^2 + (WOX5(Tf==32)-71.16)^2 +...
            (WOX5(Tf==40)-233.79)^2 + (WOX5(2880)-152.02)^2;
        % R2: SHR time course data
        R2 = (SHR(Tf==8)-186.85)^2 + (SHR(Tf==16)-174.41)^2 +...
            (SHR(Tf==24)-127.87)^2 + (SHR(Tf==32)-148.53)^2 +...
            (SHR(Tf==40)-96.60)^2 + (SHR(2880)-70.81)^2;
        % R3: AN3 time course data
        R3_cei = (AN3_cei(Tf==8)-293.69)^2 + (AN3_cei(Tf==16)-179.66)^2 +...
            (AN3_cei(Tf==24)-195.13)^2 + (AN3_cei(Tf==32)-302.27)^2 +...
            (AN3_cei(Tf==40)-561.32)^2 + (AN3_cei(2880)-228.17)^2;
        R3_qc = (AN3_qc(Tf==8)-129.33)^2 + (AN3_qc(Tf==16)-79.12)^2 +...
            (AN3_qc(Tf==24)-85.93)^2 + (AN3_qc(Tf==32)-133.11)^2 +...
            (AN3_qc(Tf==40)-247.19)^2 + (AN3_qc(2880)-100.48)^2;
        R3_endo = (AN3_endo(Tf==8)-83.86)^2 + (AN3_endo(Tf==16)-51.30)^2 +...
            (AN3_endo(Tf==24)-55.72)^2 + (AN3_endo(Tf==32)-86.31)^2 +...
            (AN3_endo(Tf==40)-160.27)^2 + (AN3_endo(2880)-65.15)^2;
        R3 = R3_cei + R3_qc + R3_endo;
        % R4: SCR time course data
        R4_endo = (SCR_endo(Tf==8)-143.88)^2 + (SCR_endo(Tf==16)-154.58)^2 +...
            (SCR_endo(Tf==24)-138.73)^2 + (SCR_endo(Tf==32)-143.50)^2 +...
            (SCR_endo(Tf==40)-150.26)^2 + (SCR_endo(2880)-84.57)^2;
        R4_cei = (SCR_cei(Tf==8)-51.78)^2 + (SCR_cei(Tf==16)-55.63)^2 +...
            (SCR_cei(Tf==24)-49.93)^2 + (SCR_cei(Tf==32)-51.63)^2 +...
            (SCR_cei(Tf==40)-54.08)^2 + (SCR_cei(2880)-30.44)^2;
        R4 = R4_endo + R4_cei;
        % R5: CYCD6 time course data
        R5 = (CYCD6(Tf==8)-9.98)^2 + (CYCD6(Tf==16)-8.85)^2 +...
            (CYCD6(Tf==24)-8.68)^2 + (CYCD6(Tf==32)-15.12)^2 +...
            (CYCD6(Tf==40)-3.25)^2 + (CYCD6(2880)-9.34)^2;
        % Get the sum of squared residuals
        R = R4 + R5;
    end

end
% ODEs to run parameter estimation

function dy = gif1_dy_est(t, y, params)
%% Initialization
dy = zeros(size(y, 1), 1);

% Initial conditions
% QC
WOX5_qc = y(1); SHR_qc = y(2); AN3_qc = y(3); SCR_qc = y(4); SSC_qc = y(5);

% CEI
SHR_cei = y(7); AN3_cei = y(8); SCR_cei = y(9);
SSC_cei = y(10); CYCD6_cei = y(11);

% VASC
WOX5_vasc = y(12); SHR_vasc = y(13);

% ENDO
SHR_endo = y(14); AN3_endo = y(15); SCR_endo = y(16); SSC_endo = y(17);

%% Parameters
% Estimated
k_3_endo = params(1);
d_3_endo = params(2);
d_5_cei = params(3);
d_8_endo = params(4);
d_1_vasc = params(5);
d_4_endo = params(6);
d_4_qc = params(7);
k_4_vasc = params(8);
d_4_vasc = params(9);
d_4_cei = params(10);
d_8_cei = params(11);

% Fitted to time course data
k_1_qc = params(12);
k_2_cei = params(13); k_2_qc = k_2_cei; k_2_endo = k_2_cei;

% Set
k_3_qc = 268.57;
k_8_qc = 5;
K_D2_qc = 500; K_D3_qc = 600; K_D4_qc = 1000;
d_3_qc = 1; d_8_qc = 0.5;
k_3_cei = 77.86; k_5_cei = 25000; k_8_cei = 5;
K_D2_cei = 500; K_D3_cei = 600; K_D4_cei = 1000;
d_3_cei = 0.5;
K_D1_vasc = 3;
k_8_endo = 5;
K_D2_endo = 1000; K_D3_endo = 600; K_D4_endo = 1000;

% Known
a_qc = 3.3099;
a_cei = 4.6875;
b_qc = 14.6346;

%% Differential equations
% QC
dy(1) = (k_1_qc * WOX5_qc) - (b_qc * WOX5_qc); % WOX5_qc
dy(2) = (a_qc * SHR_vasc) - (d_4_qc * SHR_qc); % SHR_qc
dy(3) = k_2_qc * AN3_qc; % AN3_qc
dy(4) = k_3_qc * (((K_D4_qc * SCR_qc + SSC_qc) /...
    (K_D3_qc * K_D4_qc + K_D3_qc * SHR_qc +...
    K_D4_qc * SCR_qc + SSC_qc)) + (AN3_qc /...
    (K_D2_qc + AN3_qc))) - (d_3_qc * SCR_qc); % SCR_qc
dy(5) = (k_8_qc * SHR_qc * SCR_qc) - (d_8_qc * SSC_qc); % SSC_qc

% CEI
dy(6) = 0; % WOX5_cei removed from model
dy(7) = (a_cei * SHR_vasc) - (d_4_cei * SHR_cei); % SHR_cei
dy(8) = k_2_cei * AN3_cei; % AN3_cei
dy(9) = k_3_cei * (((K_D4_cei * SCR_cei + SSC_cei) /...
    (K_D3_cei * K_D4_cei + K_D3_cei * SHR_cei +...
    K_D4_cei * SCR_cei + SSC_cei)) + (AN3_cei /...
    (K_D2_cei + AN3_cei))) - (d_3_cei * SCR_cei); % SCR_cei
dy(10) = (k_8_cei * SHR_cei * SCR_cei) - (d_8_cei * SSC_cei); % SSC_cei
dy(11) = k_5_cei * (SSC_cei / (K_D3_cei * K_D4_cei...
    + K_D4_cei * SCR_cei + K_D3_cei * SHR_cei...
    + SSC_cei)) - (d_5_cei * CYCD6_cei); % CYCD6_cei

% VASC
dy(12) = b_qc * WOX5_qc - d_1_vasc * WOX5_vasc; % WOX5_vasc
dy(13) = k_4_vasc * (K_D1_vasc / (K_D1_vasc + WOX5_vasc)) -...
    (d_4_vasc * SHR_vasc) - (a_qc * SHR_vasc) -...
    (2 * a_cei * SHR_vasc); % SHR_vasc

% ENDO
dy(14) = (a_cei * SHR_vasc) - (d_4_endo * SHR_endo); % SHR_endo
dy(15) = k_2_endo * AN3_endo; % AN3_endo
dy(16) = k_3_endo * (((K_D4_endo * SCR_endo + SSC_endo) /...
    (K_D3_endo * K_D4_endo + K_D3_endo * SHR_endo +...
    K_D4_endo * SCR_endo + SSC_endo)) + (AN3_endo /...
    (K_D2_endo + AN3_endo))) - (d_3_endo * SCR_endo); % SCR_endo
dy(17) = (k_8_endo * SHR_endo * SCR_endo) -...
    (d_8_endo * SSC_endo); % SSC_endo

end
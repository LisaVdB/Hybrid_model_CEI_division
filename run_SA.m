% Run sensitivity analysis i number of times and average results

%% Initialization
S_results = cell(1, 10);
ST_results = cell(1, 10);

%% Run SA
for i = 1:10
    [S, ST] = gif1_sensitivity_analysis(150, 17);
    S_results{i} = S;
    ST_results{i} = ST;
end

%% Combine and average results
S_results_stacked = cat(3, S_results{:});
ST_results_stacked = cat(3, ST_results{:});

S_average = mean(S_results_stacked, 3);
ST_average = mean(ST_results_stacked, 3);

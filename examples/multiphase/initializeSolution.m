function sol = initializeSolution(data)
    
    sol = struct;
    % Initial guesses for state, control, and time grid
    sol.time = initialTimeGrid(data);
    sol.states = initialStateGuess(sol,data);
    sol.controls = initialControlGuess(sol,data);
    sol.states = initialStateGuess_DF(sol,data);
    sol.params = initialParamGuess(sol,data);
    % load("data_1.mat","solution","xsim");
    % sol.time = solution.time;
    % sol.states = xsim;
    % sol.controls = solution.controls;
    % sol.params = solution.params;
end
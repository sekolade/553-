root = 'C:\Users\savka\Desktop\ders\5.1\553\notes';
addpath(genpath(root));



% {disp_iter, c, g_n, ss, mu, stop_val, total_iter}


% % Objective
% f = @(y,par) (-2*y(1) + y(2));
% 
% % h(y)=0
% h_eq = {@(y,par) (y(2) - y(1)^2)};
% g_ineq = {};
% 
% prob_opts.mode = 'penalty';
% prob_opts.W = 1;   % SPD scalar
% prob = build_problem(f, g_ineq, h_eq, prob_opts);
% 
% p = struct();
% p.prob = prob;
% p.params = struct();
% 
% p.method = 'Newton';
% p.max_iter = 100;
% p.tol = 1e-6;
% p.stop_type = 5;     % ||grad||
% p.mu_min = 1; p.mu_max = 100;   % penalty büyüsün istersen
% p.delta = 1e-6;
% 
% p.iiu_c = 1;
% p.linesearch = 'Parabolic';
% p.ls_tol = 1e-4;
% p.newton_1d_iter = 20;
% 
%%




% y0 = [1;1]; MUST be strictly feasible for barrier
% [y_sol, ch, t_data] = run_optimization_routine(y0, p)
% 



%%






% --- Q5: min x1 - x2 + x2^2  s.t. x1>=0, x2>=0
f = @(y,u) (y(1) - y(2) + y(2)^2);

% Convert to g(y)<=0:  -x1<=0, -x2<=0
g_ineq = {@(y,u) (-y(1)); @(y,u) (-y(2))};
h_eq = {};

opts = struct('mode','barrier','barrier_type','log');
prob = build_problem(f, g_ineq, h_eq, opts);

p = struct();
p.prob = prob;
p.user = [];

p.method = 'DFP';
p.max_iter = 2500;
p.tol = 1e-5;
p.stop_type = 2;

% Decrease mu toward 0 (linear schedule allows mu_min > mu_max)
p.mu_min = 1;
p.mu_max = 1e-4;

p.linesearch = 'Dichtomous';
p.ls_tol = 1e-4;
p.iiu_c = 1;
p.delta = 1e-1;
p.newton_1d_iter = 2;

p.eps_grad = 1e-6;
p.eps_hess = 1e-5;
%%

y0 = [4;4];   % MUST be strictly feasible for barrier


%%
[y_sol, ch, tdata] = run_optimization_routine(y0, p)












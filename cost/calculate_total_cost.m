function J = calculate_total_cost(y, mu, params)
%CALCULATE_TOTAL_COST Total cost evaluator:
%   J(y) = f(y) + mu * Penalty(y)           (penalty)
%   J(y) = f(y) + mu * Barrier(y)           (barrier)
%   J(y) = f(y) + mu*Penalty(y) + mu*Barrier(y) (mixed)
%   J(y) = f(y)                             (none)
%
% Required in params:
%   params.prob   : created by build_problem(...)
%   params.user   : (optional) arbitrary user parameters for f/g/h handles
%
% Notes:
% - Constraint handles are always called as handle(y, params.user).

    if ~isfield(params,'prob') || isempty(params.prob)
        error('params.prob must be provided (use build_problem).');
    end

    prob = params.prob;
    user = [];
    if isfield(params,'user'), user = params.user; end

    % Base objective
    f = prob.f(y, user);

    switch lower(prob.opts.mode)
        case 'none'
            J = f;

        case 'penalty'
            P = penalty_eval(y, prob, user);
            J = f + mu * P;

        case 'barrier'
            B = barrier_eval(y, prob, user);
            J = f + mu * B;

        case 'mixed'
            P = penalty_eval(y, prob, user);
            B = barrier_eval(y, prob, user);
            J = f + mu * P + mu * B;

        otherwise
            error('Unknown prob.opts.mode.');
    end
end

% =========================
% Penalty / Barrier helpers
% =========================

function P = penalty_eval(y, prob, user)
% Quadratic penalty with SPD matrices:
%   P = 1/2 * g_plus' * Gamma * g_plus + 1/2 * h' * W * h
% where g_plus = max(g,0) penalizes only violations.

    g = eval_g_ineq(y, prob, user); % p x 1
    h = eval_h_eq(y, prob, user);   % q x 1

    if isempty(g)
        Pg = 0;
    else
        gplus = max(g, 0);
        Pg = 0.5 * (gplus.' * prob.opts.Gamma * gplus);
    end

    if isempty(h)
        Ph = 0;
    else
        Ph = 0.5 * (h.' * prob.opts.W * h);
    end

    P = Pg + Ph;
end

function B = barrier_eval(y, prob, user)
% Simple interior-point barrier:
% g_i(y) < 0 required. Otherwise B = +Inf.

    g = eval_g_ineq(y, prob, user);
    if isempty(g)
        B = 0;
        return;
    end

    if any(g > 0)
        B = Inf;                 % outside the domain
        return;
    end

    switch lower(prob.opts.barrier_type)
        case 'log'
            B = -sum(log(-g));    % classic log barrier
        case 'inv'
            B = sum(1./(-g));     % inverse barrier (less common than log)
        otherwise
            error('Unknown barrier_type.');
    end
end

function g = eval_g_ineq(y, prob, user)
    p = numel(prob.g_ineq);
    if p == 0, g = []; return; end
    g = zeros(p,1);
    for i=1:p
        g(i) = prob.g_ineq{i}(y, user);
    end
end

function h = eval_h_eq(y, prob, user)
    q = numel(prob.h_eq);
    if q == 0, h = []; return; end
    h = zeros(q,1);
    for j=1:q
        h(j) = prob.h_eq{j}(y, user);
    end
end

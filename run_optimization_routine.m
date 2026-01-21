function [y, ch, t_data] = run_optimization_routine(y, p)
%RUN_OPTIMIZATION_ROUTINE Core optimization loop (no GUI).
%
% Required algorithm fields in p:
%   p.method    : 'SD' | 'Newton' | 'FR' | 'DFP'
%   p.max_iter  : maximum iterations
%   p.tol       : stopping tolerance
%   p.stop_type : 1|2|3|4|5 (same meaning as your original)
%   p.mu_min, p.mu_max : mu schedule endpoints
%
% Line-search / Newton options in p:
%   p.linesearch      : 'Fibonacci'|'Dichotomous'|'Parabolic'|'Newton1D'
%   p.ls_tol          : line-search tolerance
%   p.iiu_c           : initial interval constant
%   p.delta           : minimum eigenvalue threshold for modified Newton
%   p.newton_1d_iter  : max iter for 1D Newton line search
%
% Problem definition in p:
%   p.prob (from build_problem), p.user (optional), and eps_grad/eps_hess (optional)

    ch = [];
    n = length(y);

    disp_iter  = 0;   % logging counter
    total_iter = 0;   % actual iterations

    % State
    st.n = n;
    st.iter = 0;

    % DFP state
    if strcmp(p.method,'DFP')
        st.S = eye(n);
        st.j = 1;
    end

    % FR state
    if strcmp(p.method,'FR')
        st.g_old = [];
        st.d_old = [];
        st.j = 1;
    end

    data_buffer = {};
    c_prev = calculate_total_cost(y, p.mu_max, p);

    for k = 1:p.max_iter

        % Mu schedule (linear, same structure as your original code)
        if p.max_iter > 1
            mu = p.mu_min + (p.mu_max - p.mu_min) * (k-1) / (p.max_iter - 1);
        else
            mu = p.mu_max;
        end

        % Execute one step
        log_now = false;
        if strcmp(p.method,'SD')
            [y2, st, c, g_n, ss] = execute_sd(y, st, mu, p);
            log_now = true;
        elseif strcmp(p.method,'Newton')
            [y2, st, c, g_n, ss] = execute_newton(y, st, mu, p);
            log_now = true;
        elseif strcmp(p.method,'FR')
            [y2, st, c, g_n, ss, log_now] = execute_fr(y, st, mu, p);
        elseif strcmp(p.method,'DFP')
            [y2, st, c, g_n, ss, log_now] = execute_dfp(y, st, mu, p);
        else
            error('Unknown method.');
        end

        total_iter = total_iter + 1;

        % Stopping value (same semantics)
        if p.stop_type == 1
            stop_val = abs(c - c_prev);
        elseif p.stop_type == 2
            stop_val = ss;
        elseif p.stop_type == 3
            stop_val = abs(c - c_prev) / (abs(c_prev) + 1e-10);
        elseif p.stop_type == 4
            stop_val = ss / (norm(y) + 1e-10);
        elseif p.stop_type == 5
            stop_val = g_n;
        else
            error('Unknown stop_type.');
        end

        c_prev = c;
        is_converged = (stop_val < p.tol);

        % Logging rule: first step OR reset-trigger OR converged
        if log_now || k == 1 || is_converged
            disp_iter = disp_iter + 1;

            ch = [ch; c];
            data_buffer(disp_iter,:) = {disp_iter, c, g_n, ss, mu, stop_val, total_iter};

            if is_converged
                y = y2;
                t_data = data_buffer;
                return;
            end
        end

        y = y2;
    end

    t_data = data_buffer;
end

function [y_new, state, cost, g_norm, step_size] = execute_newton(y, state, mu, params)
%EXECUTE_NEWTON Modified Newton step:
% - Computes numerical Hessian
% - Shifts it to be positive definite if min eigenvalue < delta
% - Performs line search along d = -H^{-1} g

    g = calculate_gradient(y, mu, params);
    cost = calculate_total_cost(y, mu, params);
    g_norm = norm(g);

    H = calculate_hessian_wrapper(y, mu, params);

    evals = eig(H);
    min_eig = min(evals);
    delta_val = params.delta;

    if min_eig < delta_val
        epsilon = delta_val - min_eig;
        H = H + epsilon * eye(state.n);
    end

    d = -H \ g;

    alpha = perform_line_search(y, d, mu, params);
    y_new = y + alpha*d;
    step_size = norm(y_new - y);

    state.iter = state.iter + 1;
end

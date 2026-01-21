function [y_new, state, cost, g_norm, step_size] = execute_sd(y, state, mu, params)
%EXECUTE_SD Steepest descent step with line search.

    g = calculate_gradient(y, mu, params);
    cost = calculate_total_cost(y, mu, params);
    g_norm = norm(g);

    d = -g;
    alpha = perform_line_search(y, d, mu, params);

    y_new = y + alpha*d;
    step_size = norm(y_new - y);

    state.iter = state.iter + 1;
end

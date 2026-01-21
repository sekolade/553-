function alpha_star = perform_line_search(y, d, mu, params)
%PERFORM_LINE_SEARCH Choose a line-search method.

    [a, b, t_mid] = initial_interval_of_uncertinity(y, d, mu, params);

    switch params.linesearch
        case 'Fibonacci'
            alpha_star = fibonacci_search(y, d, mu, params, a, b);
        case 'Dichotomous'
            alpha_star = dichotomous_search(y, d, mu, params, a, b);
        case 'Parabolic'
            alpha_star = parabolic_fit(y, d, mu, params, a, b, t_mid);
        case 'Newton1D'
            alpha_star = newton_1d_search(y, d, mu, params, t_mid);
        otherwise
            alpha_star = (a + b) / 2;
    end
end

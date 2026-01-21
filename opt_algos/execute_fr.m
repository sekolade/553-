function [y2, st, c, gn, ss, log_now] = execute_fr(y, st, mu, params)
%EXECUTE_FR Fletcher-Reeves conjugate gradient with periodic reset.

    g = calculate_gradient(y, mu, params);
    c = calculate_total_cost(y, mu, params);
    gn = norm(g);

    log_now = false;

    if st.iter == 0
        d = -g;
        st.j = 1;
        log_now = true;  % log the first step
    else
        limit = st.n;

        if st.j < limit
            g_prev = st.g_old;
            beta = (g' * g) / (g_prev' * g_prev);
            d = -g + beta * st.d_old;
            st.j = st.j + 1;
        else
            % Reset
            d = -g;
            st.j = 1;
            log_now = true; % log reset step
        end
    end

    st.g_old = g;
    st.d_old = d;
    st.iter = st.iter + 1;

    a = perform_line_search(y, d, mu, params);
    y2 = y + a*d;
    ss = norm(y2 - y);
end

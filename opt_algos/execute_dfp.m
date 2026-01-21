function [y2, st, c, gn, ss, log_now] = execute_dfp(y, st, mu, params)
%EXECUTE_DFP DFP quasi-Newton (inverse Hessian approximation) with periodic reset.

    g = calculate_gradient(y, mu, params);
    c = calculate_total_cost(y, mu, params);
    gn = norm(g);

    n = length(y);
    if st.iter == 0
        st.S = eye(n);
        st.j = 1;
    end

    d = -st.S * g;

    a = perform_line_search(y, d, mu, params);
    y2 = y + a*d;
    ss = norm(y2 - y);

    limit = n/2 + 1;
    log_now = false;

    if st.j < limit
        s = y2 - y;
        g_next = calculate_gradient(y2, mu, params);
        q = g_next - g;

        % DFP update with a safeguard
        if (s' * q) > 1e-10
            term1 = (s * s') / (s' * q);
            term2 = (st.S * q * q' * st.S) / (q' * st.S * q);
            st.S = st.S + term1 - term2;
        end

        st.j = st.j + 1;
    else
        % Reset
        st.S = eye(n);
        st.j = 1;
        log_now = true;
    end

    st.iter = st.iter + 1;
end

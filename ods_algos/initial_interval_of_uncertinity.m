function [a, b, t2] = initial_interval_of_uncertinity(y, d, mu, params)
%INITIAL_INTERVAL_OF_UNCERTINITY Bracket a minimum along direction d.

    g_curr = calculate_gradient(y, mu, params);
    g_mag  = norm(g_curr);

    alpha_max = params.iiu_c / g_mag;
    if isfield(params,'lambda_max_current') && alpha_max > params.lambda_max_current
        alpha_max = params.lambda_max_current;
    end

    t1 = 0;
    t3 = alpha_max;
    t2 = t3/2;

    f1 = cost_on_line(t1, y, d, mu, params);
    f2 = cost_on_line(t2, y, d, mu, params);
    f3 = cost_on_line(t3, y, d, mu, params);

    % Shrink interval if not decreasing at mid
    %while f2 >= f1 || f3<1e5
    while f2 >= f1
        alpha_max = alpha_max / 2;
        t1 = 0; t3 = alpha_max; t2 = t3/2;

        f1 = cost_on_line(t1, y, d, mu, params);
        f2 = cost_on_line(t2, y, d, mu, params);
        f3 = cost_on_line(t3, y, d, mu, params);
        
        if alpha_max < 1e-8, break; end
    end

    % Expand interval if still decreasing towards t3
    while f2 >= f3 
        if isfield(params,'lambda_max_current') && alpha_max > params.lambda_max_current
            t3= params.lambda_max_current;
            t2=t3/2;
            break;
        end
        %if f3>1e5, break; end

        alpha_max = alpha_max * 2;
        t1 = 0; t3 = alpha_max; t2 = t3/2;


        f1 = cost_on_line(t1, y, d, mu, params);
        f2 = cost_on_line(t2, y, d, mu, params);
        f3 = cost_on_line(t3, y, d, mu, params);

        if alpha_max > 1e5, break; end


    end

    a = t1;
    b = t3;
end

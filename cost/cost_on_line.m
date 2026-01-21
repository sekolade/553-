function val = cost_on_line(alpha, y, d, mu, params)
%COST_ON_LINE Evaluate phi(alpha) = J(y + alpha d).





    
    val = calculate_total_cost(y + alpha*d, mu, params);



end

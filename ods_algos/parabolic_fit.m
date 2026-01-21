function alpha = parabolic_fit(y, d, mu, params, a, b, t_mid)
    t1 = a; t2 = t_mid; t3 = b;
    f1 = cost_on_line(t1, y, d, mu, params);
    f2 = cost_on_line(t2, y, d, mu, params);
    f3 = cost_on_line(t3, y, d, mu, params);
    
    M = [1, t1, t1^2; 1, t2, t2^2; 1, t3, t3^2];
    rhs = [f1; f2; f3];
    coeffs = M \ rhs;
    c1 = coeffs(2); c2 = coeffs(3);
    
    if abs(c2) > 1e-10
        alpha = -c1 / (2 * c2);
    else
        alpha = t2; 
    end
    if alpha < a || alpha > b
        alpha = t2; end
end
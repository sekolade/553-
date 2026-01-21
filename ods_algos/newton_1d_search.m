function alpha = newton_1d_search(y, d, mu, params,b)
    alpha = b;
    tol = params.ls_tol; 
    max_iter = params.newton_1d_iter; 
    
    for i = 1:max_iter
        y_curr = y + alpha * d;
        g = calculate_gradient(y_curr, mu, params);
        H = calculate_hessian_wrapper(y_curr, mu, params); 
        
        phi_prime = g' * d;
        phi_double_prime = d' * H * d;
        
        if abs(phi_double_prime) < 1e-8, break; end
        
        alpha_new = alpha - phi_prime / phi_double_prime;
        
        if abs(alpha_new - alpha)  < tol 
            alpha = alpha_new; 
            break; 
        end
        alpha = alpha_new;
    end
    if alpha < 0, alpha = 1e-5; end
end
function H = calculate_hessian_wrapper(y, mu, params)
%CALCULATE_HESSIAN_WRAPPER Returns numerical Hessian of total cost.
    H = numerical_hessian(y, mu, params);
end

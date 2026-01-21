function g = numerical_gradient(y, mu, params)
%NUMERICAL_GRADIENT Central difference numerical gradient:
%   g(i) = (J(y + h e_i) - J(y - h e_i)) / (2h)
%
% Optional:
%   params.eps_grad : base FD step (default 1e-6)

    if isfield(params,'eps_grad') && ~isempty(params.eps_grad)
        eps_step = params.eps_grad;
    else
        eps_step = 1e-6;
    end

    n = length(y);
    g = zeros(n,1);

    for i=1:n
        h = eps_step * max(1, abs(y(i)));

        y_plus  = y; y_plus(i)  = y_plus(i)  + h;
        y_minus = y; y_minus(i) = y_minus(i) - h;

        g(i) = (calculate_total_cost(y_plus,  mu, params) - ...
                calculate_total_cost(y_minus, mu, params)) / (2*h);
    end
end

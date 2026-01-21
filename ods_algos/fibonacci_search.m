function val = fibonacci_search(y, d, mu, p, a, b)
    eps = p.ls_tol;
    len = b - a;
    
    F = [1 1];
    while (len / F(end)) > eps
        F(end+1) = F(end) + F(end-1);
    end
    
    n = length(F); 
    if n < 3
        val = (a + b) / 2;
        return;
    end
    k = n; 
    
    x1 = a + (F(k-2)/F(k)) * (b-a);
    x2 = a + (F(k-1)/F(k)) * (b-a);
    
    c1 = cost_on_line(x1, y, d, mu, p);
    c2 = cost_on_line(x2, y, d, mu, p);
    
    while k > 3
        k = k - 1;
        dist = b - a;
        
        if c1 < c2
            b = x2;
            x2 = x1;
            c2 = c1;
            x1 = a + (F(k-2)/F(k)) * dist;
            c1 = cost_on_line(x1, y, d, mu, p);
        else
            a = x1;
            x1 = x2;
            c1 = c2;
            x2 = a + (F(k-1)/F(k)) * dist;
            c2 = cost_on_line(x2, y, d, mu, p);
        end
    end
    val = (a + b) / 2;
end
function val = dichotomous_search(y, d, mu, p, a, b)
    tol = p.ls_tol; 
    e = tol / 10;
   
    while (b - a) >= tol
        mid = (a + b) / 2;
        
        x1 = mid - e;
        x2 = mid + e;
        
        f1 = cost_on_line(x1, y, d, mu, p);
        f2 = cost_on_line(x2, y, d, mu, p);
        
        if f1 < f2
            b = x2;
        else
            a = x1;
        end
    end
    
    val = (a + b) / 2;
end


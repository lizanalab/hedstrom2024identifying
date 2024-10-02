function val = hdm_ana(d, n, a)
    if d == 0
        %val = (a^n);
        val = 0; % inspired by regular power-law dropoff where distance is zero and contacts are thus ill-defined
        return
    end
    if d >= 2^n
        val = hdm_ana(2^n-1, n, a);
        return
    end
    g = @(x) (x <= 0.5) .* ((x < 0) .* 0 + (x >= 0) .* x) + (x > 0.5) .* ((x > 1) .* 0 + (x <= 1) .* (1-x));
    b = 1/2-a;
    val = 0.0;
    for k = 0:n-1
        val = val + 2^n*g(d*2.^(k-n))*a^k*b.^(1 - eq(k, n))/(4^max(n-k-1, 0));
    end
    val = val/(2^n-d);
end
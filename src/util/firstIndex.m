function idx = firstIndex(a, value)
    tol = 1e-6;
    idx = find(abs(a - value) < tol, 1, 'first');
end
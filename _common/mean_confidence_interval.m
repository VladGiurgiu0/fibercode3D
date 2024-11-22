function ci = mean_confidence_interval(y,alpha) % m,s,n,d are local variables
    n = numel(y);
    m = mean(y,'omitnan');
    s = std(y,'omitnan') / sqrt(n);
    d = s * -tinv(alpha/2, max(0,n-1));
    ci = [m-d, m+d];
end
function S_bar = pf_weight(S_bar,z,params)
    S_bar.W = sensor_model(z,params.Sigma_Q,S_bar.X(1:2,:),params.thresh_avg_likelihood);
end

function p = sensor_model(z, Sigma_Q, X,thresh_avg_likelihood)

p = exp(-0.5 * ((z(1) - X(1,:)).^2/Sigma_Q(1) + (z(2) - X(2,:)).^2 / Sigma_Q(4)));

if mean(p) < thresh_avg_likelihood %probably an outlier, or bad parameter settings
    p(:) = 1;
end
p = p / sum(p);
end
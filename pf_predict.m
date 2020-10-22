function S_bar = pf_predict(S,params)
    N = size(S.X,1);
    M = size(S.X,2);
    
    S_bar.X = S.X;
    
    S_bar.X= S_bar.X + randn(N,params.M) .* repmat(sqrt(diag(params.Sigma_R)),1,M); % Diffusion, assuming an uncorrelated sigma_R
    S_bar.W = S.W;
end
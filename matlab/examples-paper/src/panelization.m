function [t, w] = panelization(npan, nquad)
    [X, W] = legendre.gauss(nquad); 
    X = (X+1)/2; W = W/2;   % on [0,1]
    t = bsxfun(@plus, X/npan, (0:npan-1)/npan);
    t = t(:);             % G-L panelization of param    
    w = repmat(W/npan, 1, npan);
    w = w(:);
end

%{
% distSample - generate a random sample from an arbitrary discrete/finite distribution
%
% how is it different from distGen? well it makes use of matlab randn function. 
%
% x - distribution function domain
% p - probability vector associated with x
% n - number of samples to be generated
% rs - samples
%
% example:

rs = distSample([1 5 8], [0.1 0.6 0.3], 1000);
hist(rs);

%
% Dan I. Malta 2009
%}
function rs = distSample(x, p, n)
    % sort
    [xs,ind]=sort(x);
    
    % normalized probability vector
    p =p / sum(p);
    
    % sort probability vector according to x
    ps = p(ind);
    
    % calculate CDF
    F = cumsum(ps);
    
    % create samples
    for i=1:n
        u = rand(1,1);
        if u<= F(1)
            rs(i) = x(1);
        elseif u > F(end-1)
            rs(i) = x(end);
        else
            ind = find(u <= F);
            rs(i) = xs(ind(1));
        end
    end
end

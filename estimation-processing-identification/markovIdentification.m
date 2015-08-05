%{
% markovIdentification utilizes a markov chain monte-carlo method to identify
% static model parameters (the model is known / given - the parameters do not).
%
% markovIdentification assumes that:
% 1) noise is additive and gaussian.
% 2) the observationg are comprised of a single and constant prior probability distribution.
%
% TODO:
% current algorithm run for 6000 iterations (a value I so was good for very
% complex models), but this should be changed. I should add a convergence
% criteria based on gelman-rubin method.
%
% y - the data 
% x0 - initial parameters (initial guess of model parameters)
% func - a function holding the suggested model
% lb - lower bounds on model parameters (if not known, choose: -inf * ones(size(x0)))
% ub - upper bounds on model parameters (if not known, choose:  inf * ones(size(x0))) 
%
% example:

%define a forward model (here y=a*exp(-bx))
myfun=@(x,c) (exp(-x(1)*c)+x(2));

% generate some noisy data
true_x = [1;2];
c=linspace(1,10,100);
y=myfun(true_x,c) + .05*randn(1,100);

% estimate parameters
x0=[5; 10];

model = markovIdentification(y,x0,@(x)(myfun(x,c)), -100 * ones(size(x0)), 100 * ones(size(x0)));

clc;
disp(['model is:                          ', num2str(true_x(1)), ' + exp(', num2str(true_x(2)), ' * C)']);
disp(['estimated model is: ', num2str(mean(model(:,1))), ' + exp(', num2str(mean(model(:,2))), ' * C)']);

%
% Dan I. Malta 2014
%}
function model = markovIdentification(y, x0, func, lb, ub)
  % housekeeping
  N = length(y);
  njumps = 5000;
  burnin = 1000;
  update = 20;
  sampleevery = 10;
  
  % initial step
  p = x0;
  s = func(x0);
  e = (N / 2) * log(norm(y - s));
  acc = zeros(1, length(p));
  rej = zeros(1, length(p));
  keepp = zeros(njumps, length(p));
  prop = ones(1, length(p));
  
  % markov chain monte carlo iterations
  for i = 1 : njumps + burnin
    % markov
    for k = 1 : length(p)
      oldp = p;
      p(k) = p(k) + randn * prop(k);
      
      % bounds
      if p(k) < lb(k) || p(k) > ub(k)
        p(k) = oldp(k);
        rej(k) = rej(k) + 1;
      else
        s = func(p);
        olde = e;
        e = (N / 2) * log(norm(y - s));
        if exp(olde - e) > rand
          acc(k) = acc(k) + 1;
        else
          p(k) = oldp(k);
          rej(k) = rej(k) + 1;
          e = olde;
        end
      end
    end
    
    keepp(i, :) = p;
    if rem(i, update) == 0
      prop = prop .* sqrt((1 + acc) ./ (1 + rej));
      acc = 0 * acc;
      rej = 0 * rej;
    end
  end  

  % output
  model = keepp(burnin + 1 : sampleevery : end, :);
end

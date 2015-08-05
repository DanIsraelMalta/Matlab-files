%{
% distGen - generate random numbers according to a defined discrete probability distribution
%
% P           - positive vector whos values form a PDF for the indices of P
% N           - output number of rows
% M           - output number of columns
% 'plot' - displays the input distribution and generated points (optional input)
% T           - NxM generated random numbers according to supplied PDF
%
% conceptual example:
% if P = [0.2, 0.4, 0.4] (note sum(P) = 1), then T can only take values of 1, 2 or 3.
% if we call distGen(P, 1, 10), then on average, the output T should
% contain two 1's, four 2's and four 3's.
%
% example:

X = -3:0.1:3;
P = 1+sin(X)+cos(X/10)+(X/50).^2;
T = distGen(P,1,1000,'plot');
Xrand = X(T);

%
% Dan I. Malta 2009
%}
function T = distGen(P, N, M, varargin)
    % normalize P
    Pnorm=[0 P] / sum(P);
    
    % create cumlative distribution
    Pcum = cumsum(Pnorm);
    
    % create random matrix
    N = round(N);
    M = round(M);
    R = rand(N, M);
    
    % output
    T = zeros(N, M);
    for i=1 : length(P)
        T(and(R > Pcum(i), R <= Pcum(i+1))) = i;
    end
    
    %if desired, output plot
    if nargin == 4
        if strcmp(varargin{1},'plot')
            Pfreq = N * M * P / sum(P);
            figure;
            hold on
            hist(T(T > 0), 1 : length(P))
            plot(Pfreq, 'r-o')
            legend('generated probability', 'input distribution');
            ylabel('Frequency')
            xlabel('P-vector Index')
            axis tight
            box on
        end
    end
end


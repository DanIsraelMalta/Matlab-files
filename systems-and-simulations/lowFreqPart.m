%{
% It is often required to express a transfer function as a sum of the term
% of the partial fraction expression a0/(s+p1) which is the dominant
% component at lower frequencies, and the rest of the terms combined into a
% single transfer function.
% This is what this function does.
%
% inum - input transfer function numerator
% iden - input transfer function denominator (maximum order is 9)
% a0 - a0 of ao/(s+p1)
% p1 - p1 of ao/(s+p1)
% onum - rest of the transfer function numerator
% oden - rest of the transfer function denominator
%
% example:

numc = [9.6082  29.262  48.975  16.016];
denc = [1  6.8631  24.466  46.749  32.488  0];
[a0, p1, onum, oden] = lowFreqPart(numc, denc);

input = tf(numc, denc);
outputLPF = tf(a0, [1, p1]);
outputRest = tf(onum, oden);
bodemag(input, outputLPF);
legend('input', ['low pass component ', num2str(a0), '/(s+', num2str(p1),')']);

%
% Dan I. Malta 2015
%}
function [a0, p1, onum, oden] = lowFreqPart(inum, iden)
    % housekeeping
    ld = length(iden);
    
    % low pass component
    [r, p, k] = residue(inum, iden);
    a0 = r(ld - 1);
    p1 = p(ld - 1);
    
    % rest of the terms
    if ld == 2
        rr = 0;
        pr = 1;
    else
        rr = [];
        pr = [];
        for i = ld : -1 : 3
            rr = [rr, r(ld - 1)]; %#ok<AGROW>
            pr = [pr, p(ld - 1)]; %#ok<AGROW>
        end
    end
    [onum, oden] = residue(rr, pr, k);
    
    % only real parts
    onum = real(onum);
    oden = real(oden);
end

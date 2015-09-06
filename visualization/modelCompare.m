%{
% given several test patterns centered root mean square distance from
% reference pattern and their correleation to this reference pattern,
% modelCompare draws a polar diagram at which the reference pattern is at
% the (0, 0), the radius is the cenetered rms distance and angle is the
% correlation.
%
% offcourse, the centered root mean square can be changed to any other sort
% of distance (but then the xlabel/ylabel must be changed).
%
% testPattern - an even  cellarray where first cell holds test pattern name as a string
%                                  and second cell holds an array whos first element is test pattern centered-RMS
%                                  and second element hold's test pattern correlation [0..1]
% h- handle to produced figure
%
% example:

h = modelCompare({'Dan', [2 , 0.99], ...
                                          'Yaron', [1.3, 0.8], ...
                                          'Oleg', [4, 0.85], ...
                                          'Ben', [4.5, 0.2]});

%
% Dan I. Malta 2015
%}
function h = modelCompare(testPattern)
    % housekeeping
    corr = [0.1 : 0.1 : 0.9, 0.95, 0.99];
    angle = 0 : 0.01 : pi/2;
    cosAngle = cos(angle);
    sinAngle = sin(angle);
    len = numel(testPattern);
    crmsMax = -Inf;
    crmsMin = Inf;
    inc = 0.1;
    
     % visualization
    h = figure('color', [1, 1, 1]);
    hold on;
    
    % plot reference & test patterns
    for i = 1 : 2 : len - 1
        Rmap = pi/2 * sqrt(1 - testPattern{i + 1}(2));
        cRmap = cos(Rmap);
        sRmap = sin(Rmap);
        plot(testPattern{i + 1}(1) * cRmap, testPattern{i + 1}(1) * sRmap, 'ob', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        text((testPattern{i + 1}(1) + inc) * cRmap, (testPattern{i + 1}(1) + inc) * sRmap, testPattern{i});
       if testPattern{i + 1}(1) > crmsMax
           crmsMax = testPattern{i + 1}(1);
       end
       if testPattern{i + 1}(1) < crmsMin
           crmsMin = testPattern{i + 1}(1);
       end
    end
    
    % plot crms-equal circles
    crmsStep = (crmsMax - crmsMin) / 10;
    crmsMax = crmsMax + crmsStep;
    crmsMin = crmsMin - crmsStep;
    crmsStep = (crmsMax - crmsMin) / 10;
    for r = crmsMin : crmsStep : crmsMax
       if abs(r - crmsMax) > 1e-6
            plot(r * cosAngle, r * sinAngle, '--k');
        else
            plot(r * cosAngle, r * sinAngle, 'k');
        end
    end
       
    % plot correleation lines
    for R = corr
        Rmap = pi/2 * sqrt(1 - R);
        line([0, crmsMax * cos(Rmap)], [0, crmsMax * sin(Rmap)], 'color', 'k', 'LineStyle', ':');
        text((crmsMax + inc) * cos(Rmap), (crmsMax + inc) * sin(Rmap), num2str(R)) ;
    end
   
    % finish
    set(gca, 'XTick', crmsMin : crmsStep : crmsMax);
    set(gca, 'YTick', crmsMin : crmsStep : crmsMax);
    xlabel('Centered Root Mean Square');
    ylabel('Centered Root Mean Square');
end

clear; clc; close all;

rho = 1.2; % [kg/m^3] air density
nu = 1.5e-5; % [m^2/s]


[allfilenames, pathname] =  uigetfile('*.*', 'MultiSelect', 'on');
allfilenames = cellstr(allfilenames);
nFiles = length(allfilenames);
fprintf('%d files to read', nFiles)
disp(allfilenames)
for i=1:nFiles
    %% Import all data 
    [~, filename, filetype] = fileparts(allfilenames{i});
    fprintf('Current file is: \n %s \n', filename)
    [xStep,yStep,trueAirspeed,attackAng,slipAng,staticP,totalP] =...
            import5HP([pathname, filename, filetype]);
   
    
    u = trueAirspeed .* cosd(attackAng) .* cosd(slipAng);
    v = trueAirspeed .* sind(slipAng);
    w = trueAirspeed .* sind(attackAng) .* cosd(slipAng);
    qP = totalP-staticP; % Dynamic Pressure
   
    xyStep = round([xStep, yStep]); % Round to the nearest step size

    % Find unique XY coordinates
    [xyStep,~, uidx] = unique(xyStep, 'rows', 'stable'); 
    xyInch = xyStep ./ 50.8625; % Convert steps to inches
    
    % Accumulate average values using XY coordinate index
    [~, data.uMean, data.uStd] = accumulateStats(u, uidx);
    [~, data.vMean, data.vStd] = accumulateStats(v, uidx);
    [~, data.wMean, data.wStd] = accumulateStats(w, uidx);
    
    [xx, yy, uF] = griddedInterp(xyInch, data.uMean, 100, 100);
    figure
    surf(xx, yy, uF(xx, yy))
    figure
    contourf(xx, yy, uF(xx,yy))
    [~, ~, vF] = griddedInterp(xyInch, data.vMean, 100, 100);
    [~, ~, wF] = griddedInterp(xyInch, data.wMean, 100, 100);
    
    
    %real did this crap
    ycol = xyInch(:,2);
    xcol = xyInch(:,1);
    counter = 1;
    for i1 = 1:length(xyInch(:,1));
        if ycol(i1) > 12
            newx(counter) = xcol(i1);
            newy(counter)= ycol(i1);
            formeanw(counter)= data.wMean(i1);
            formeanv(counter) = data.vMean(i1);
            counter = counter + 1;
        end
    end
    
    meanneww = mean(formeanw);
    meannewv = mean(formeanv);
    
    for i2 = 1:length(xyInch(:,1));
        neww(i2) = data.wMean(i2) - meanneww;
        newv(i2) = data.vMean(i2) - meannewv;
    end
    figure
    quiver(xyInch(:,1),xyInch(:,2),1*transpose(neww),1*transpose(newv))
    figure
    quiver(xyInch(:,1),xyInch(:,2),data.wMean,data.vMean)
end



function [xGroupped, xMean, xStd] = accumulateStats(x, uidx)
    nGroups = max(uidx); % Number of locations sampled
    xGroupped{nGroups, 1} = [];
    for i=1:nGroups
        gIdx = (uidx == i);
        xGroupped{i} = x(gIdx);
    end

    xMean = accumarray(uidx, x, [], @mean);
    xStd = accumarray(uidx, x, [], @std);
end


function [xx, yy, F] = griddedInterp(xy, fxy, nX, nY)
% Input:
%   xy = [x, y] coordinates
%   fxy = f(x,y)
%   nX, nY = number of x,y points desired

   
    F = scatteredInterpolant(xy(:,1),...
        xy(:,2), fxy);

    [XY_min, ~] = min(xy); % Lower X,Y limits
    X1 = XY_min(1);
    Y1 = XY_min(2);
    [XY_max, ~] = max(xy); % Upper X,Y limits
    X2 = XY_max(1);
    Y2 = XY_max(2);

    X_vec = linspace(X1, X2, nX);
    Y_Vec = linspace(Y1, Y2, nY);
    [xx, yy] = meshgrid(X_vec ,Y_Vec);
end

clear;  close all;

speed = 26;  % 20 or 26 meters/s
firstLocation = 'v'; % Spanwise loc. of first data point 'p', 'v', 'mp', 'mv' 

chord = 0.197; % [m] chord length
rho = 1.2; % [kg/m^3] air density
nu = 1.5e-5; % [m^2/s]
zRef = 0; % [Inch]
yRef = 23.5; % [Inch]

[allfilenames, pathname] =  uigetfile('*.*', 'MultiSelect', 'on');
allfilenames = cellstr(allfilenames);
nFiles = length(allfilenames);
fprintf('%d files to read', nFiles)
disp(allfilenames)
for i=1:nFiles
    %% Import all data 
    [~, filename, filetype] = fileparts(allfilenames{i});
    fprintf('Current file is: \n %s \n', filename)
    data(1) = struct(); % Initialize a struct for data
    data.filename = filename;
    data.z1 = firstLocation;
    [zStep,yStep,trueAirspeed,attackAng,slipAng,staticP,totalP] =...
            import5HP([pathname, filename, filetype]);
    % Offset the angles
%     [attackOff, slipOff ] = meanOffset(speed);
%     attackAng = attackAng - attackOff;
%     slipAng = slipAng - slipOff;
    slipAng = slipAng - mean(slipAng);
    attackAng = attackAng - mean(attackAng);
    
    u = trueAirspeed .* cosd(attackAng) .* cosd(slipAng);
    v = trueAirspeed .* sind(slipAng);
    w = trueAirspeed .* sind(attackAng) .* cosd(slipAng);
    qP = totalP-staticP; % Dynamic Pressure
   
    zyStep = round([zStep, yStep]); % Round to the nearest step size

    % Find unique ZY coordinates
    [zyStep,~, uidz] = unique(zyStep, 'rows', 'stable'); 
    zyInchRelative = zyStep ./ 50.8625; % Convert steps to inches
    data.zyInch = [zRef, yRef] - zyInchRelative; 
    data.zyInch(:,1) = data.zyInch(:,1)*(-1); % Correccts sign of z coordinate
    % Get user input data
    data = getConfigData(data);
    
    % Accumulate average/stdev values using ZY coordinate index
    [~, data.u, data.uStd] = accumulateStats(u, uidz);
    [~, data.v, data.vStd] = accumulateStats(v, uidz);
    [~, data.w, data.wStd] = accumulateStats(w, uidz);
    
    data = griddify(data);
    data = scaleData(data);
    data = getVorticity(data);
    data = getDrag(data, chord);
    data = uncert(data);
    
    quickPlots(data)
    
    data.filename = filename;
    save(['.\Wake Data\' filename '.mat'], '-struct', 'data')
end

function data = griddify(data)
% This function is only designed to work for evenly spaced grids with
% constant delta z and delta y where the deltas need not be equal
    z = data.zyInch(:,1); % z VECTOR for data locations
    y = data.zyInch(:,2); % y VECTOR for data locatoins
    uZ = unique(z); % Ordered unique z points;
    uY = unique(y);
    data.dz = uZ(2) - uZ(1); % Delta z for data locations
    data.dy = uY(2) - uY(1); % Delta y for data locations
    [data.zz, data.yy] = meshgrid(uZ, uY); 
    % xq, yq are ordered (low to high), unique grid point coordinates
    % corresponding to actual data locations
    nZ = 100;
    nY = 100;
    zi = linspace(uZ(1), uZ(end), nZ); % z for ordered, interp'ed grid
    yi = linspace(uY(1), uY(end),nY);
    data.dzi = zi(2)-zi(1); % Delta z for interp'ed grid
    data.dyi = yi(2)-yi(1); % Delta y for interp'ed grid
    [data.zzi, data.yyi] = meshgrid(zi, yi);
    % zqi, yqi are ordered coords for interpolating between data collection
    % locations.
    
    % First grid all the true data
    data.uu = griddata(z, y, data.u, data.zz, data.yy);
    data.vv = griddata(z, y, data.v, data.zz, data.yy);
    data.ww = griddata(z, y, data.w, data.zz, data.yy);
        % Grid the standard deviations
    data.uuStd = griddata(z, y, data.uStd, data.zz, data.yy);
    data.vvStd = griddata(z, y, data.vStd, data.zz, data.yy);
    data.wwStd = griddata(z, y, data.wStd, data.zz, data.yy);
    
    % Now grid with interpolation
    data.uui = griddata(z, y, data.u, data.zzi, data.yyi);
    data.vvi = griddata(z, y, data.v, data.zzi, data.yyi);
    data.wwi = griddata(z, y, data.w, data.zzi, data.yyi);    
    
end

function [xGroupped, xMean, xStd] = accumulateStats(x, uidx)
% Finds the mean and standard deviation of a the data after grouping. For
% example, data may be grouped by unique coordinates (determined by uidx)
% and the output will give the statistics at each coordinate.

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
%   fxy = f(x,y) = f([xy])
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
 
    % Example:
%     [xx, yy, data.uu] = griddedInterp(data.xyInch, data.u, 100, 100);
%     figure
%     surf(xx, yy, uF(xx, yy))
%     figure
%     contourf(xx, yy, uF(xx,yy))
%     [~, ~, vF] = griddedInterp(data.xyInch, data.v, 100, 100);
%     [~, ~, wF] = griddedInterp(data.xyInch, data.w, 100, 100);
   
end

function data = getConfigData(data)
    inptFields = {'Length (inch)', 'Width (inch)', 'Configuration'};
    fieldSizes = [1 20; 1 20; 1 20];
    
    % inptFields = {'Configuration'};
    % fieldSizes = [1 20];
    inpt = inputdlg(inptFields, data.filename, fieldSizes);
    
    data.L = str2double(inpt{1});
    data.W = str2double(inpt{2});
    data.config = inpt{3};

    % config = inpt{1};
end

function data = scaleData(data)
    data = getScalingFactors(data);
    % Scale the data
    data.yys = (data.yy - data.yMin)/data.yHalf;
    minZ = min(min(data.zz));
    maxZ = max(max(data.zz));
    if isnan(data.W) % No CVG
        data.zzs = (data.zz - minZ)/(maxZ - minZ);
    else % Yes CVG, divide by width, align peaks to 0
        switch data.z1
            case 'p' % Peak
                zOff = 0;
            case 'v' % Valley
                zOff = 0.5;
            case 'mp' % Mid then peak
                zOff= -0.25;
            case 'mv' % Mid then valley
                zOff = 0.25;
        end
        data.zzs = (data.zz - data.zz(1) + zOff)/data.W;
    end
    data.uus = data.uu/data.uInf;
    data.vvs = data.vv/data.uInf;
    data.wws = data.ww/data.uInf;
        
end

function data = getScalingFactors(data)
% This function finds the scaling factors based on the mean y profile of
% the u data.

    % First find the mean profile
    uMean = mean(data.uu,2);
    y = data.yy(:,1); % Any column of yy gives the y VECTOR

    % Find minimuum u and its y location
    [uMin, yMinIdx] = min(uMean);
    yMin = y(yMinIdx);
    data.yMin = yMin;
    
    % Find the freestream
    uInfB = mean(uMean(1:5)); % Mean of bottom 5 points
    uInfT = mean( uMean(end-4: end)); % Mean of top 5 points
    
    data.uInf = mean([uInfB uInfT]); % Mean freestream
    
    % Find the half recovery points
    uMidB = (uInfB + uMin)/2; % Bottom half recovery speed
    uMidT = (uInfT + uMin)/2; % Top half recovery speed
    uFB = griddedInterpolant(y(1:yMinIdx), uMean(1:yMinIdx));
    uFT = griddedInterpolant(y(yMinIdx:end), uMean(yMinIdx:end));
    funB = @(y) uFB(y) - uMidB;
    funT = @(y) uFT(y) - uMidT;
    yHalfB = fzero(funB, yMin); % Bottom half recovery, absolute y
    yHalfT = fzero(funT, yMin); % Top half recovery
    
    data.yHalf = mean([(yMin-yHalfB) (yHalfT-yMin)]); % Mean half recovery
    % Plot the mean profile
    figure('Name', 'Mean Profile')
    plot( uMean, y)
    xlabel('Mean u profile [m/s]')
    ylabel('y [in]')
    hold on
    x1 = min(xlim);
    fill([x1; uMean; x1], [y(1); y; y(end)], 'r', 'FaceAlpha', 0.1,...
      'LineStyle', 'none')


    yL= ylim;
    y1 = yL(1); % Lower ylim
    y2 = yL(2); % Upper ylim
    plot([uInfB uInfB], [y1 yMin], 'k', 'LineWidth', 1) % Bottom uInf
    plot([uInfT uInfT], [y2 yMin], 'k', 'LineWidth', 1) % Top uInf
    plot([data.uInf data.uInf], [y1 y2], '--k', 'LineWidth', 1.2) % uInf
    plot(xlim, [yMin yMin], '--k') % Horizontal at uMin
    plot([uMin uMin], ylim, '-.k')
    plot([uMidB, uMidT], [yHalfB, yHalfT], '*r')
    
    
end

function data = getVorticity(data)
% This function calculates the vorticity of the given velocity field.
%   Omega_x = dw/dy - dv/dz [1/s]
    dz = data.dz * 0.0254; % Convert inches to meters
    dy = data.dy * 0.0254;
    [dwdz, dwdy] = gradient(data.ww, dz, dy);
    [dvdz, dvdy] = gradient(data.vv, dz, dy); 
    data.omegaX = dwdy - dvdz;
    
    data.maxOmega = max(max(data.omegaX));
    
    % Find the mean vorticity of the 98% wake region
    idx = data.uus <= 0.98;
    data.meanOmega = mean(abs(data.omegaX(idx)));
end

function data = getDrag(data, chord)
    n = size(data.zz,2); % The number of spanwise locations
    cd(1, n) = 0; % Initialize coefficent of drag
    cdUncert(1, n) = 0;
    for i = 1:n
       y = data.yy(:, i) * 0.0254; % [m] The y vector for the current spanwise location.
       % Note that all y vectors should likely be the same
       us = data.uus(:, i);
       integrand = us.*(1-us);
       cd(i) = (2/chord) * trapz(y, integrand);
       % Calc uncert for cd
       usUncert = (data.uuStd(:, i)/data.uInf).*us;
       intUncert = usUncert.*(1-usUncert);
       cdUncert(i) = (2/chord) * trapz(y, intUncert);
    end
    data.cd = cd;
    data.cdUncert = cdUncert;
end

function data = uncert(data)
     
end

function quickPlots(data)
    figure('Name', 'Contour', 'units', 'normalized',...
        'outerposition', [0 .67 .33 .33])
    contourf(data.zz, data.yy, data.uu)
    hold on 
    quiver(data.zz, data.yy, data.ww, data.vv)
    
    figure('Name', 'Scaled Contour', 'units', 'normalized',...
        'outerposition', [.33 .67 .33 .33])
    contourf(data.zzs, data.yys, data.uus)
    hold on 
    quiver(data.zzs, data.yys, data.wws, data.vvs)  
    
    figure('Name', 'Vorticity', 'units', 'normalized',...
        'outerposition', [.33 .33 .33 .33])
    contourf(data.zzs, data.yys, data.omegaX) 
    
    figure('Name', 'Vector Field', 'units', 'normalized',...
        'outerposition', [0 .33 .33 .33])
    quiver(data.zz, data.yy, data.ww, data.vv)
    
    figure('Name', 'Secction Drag Coefficient', 'units', 'normalized',...
    'outerposition', [.66 .33 .33 .33])
    errorbar(data.zzs(1,:), data.cd, data.cdUncert)
    xlabel('Scaled Span')
    ylabel('C_d')
end

function [attackOff, slipOff ]= meanOffset(speed)
    if speed == 20
        attackOff = 0.1614;
        slipOff = -0.5191;
    elseif speed == 26
        attackOff = 0.0493;
        slipOff = -0.5967;
    end       
end


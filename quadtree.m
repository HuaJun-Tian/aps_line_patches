function [nb, err, npoints, centers, values, polys, xlims, ylims] = quadtree(xy,val,thresh,maxlevel,startlevel)

% Quadtree subsampling of scattered data modified after Decriem (2009) and
% Gonzalez (2015).
%
% Usage: [nb, err, npoints, centers, values, polys, xlims, ylims] = quadtree(xy,val,thresh,maxlevel,startlevel)
% Input Parameters:
%   xy    : Mx3 array of pixels positions
%   val   : Mx1 array of the value of the pixels
%   thresh : acceptance  threshold. If variance(polygon)>threshold
%           then the polygon is subdivided.
%   maxlevel  : maximum number of subdivisions.
%   startlevel : start level (usually = 1).
%
% Output Parameters:
%   nb: number of polygons
%   err: sum of the errors
%   npoints: number of pixels contained in the quad-tree
%   centers: nx2 array of the center of mass of each polygon
%   vavlues: n   array of the value affected to each polygon
%   polys: cell{n} array of cells containing the polygons coordinates
%   xlims: nx2 array of x-coordinate limits for each selected polygon
%   ylims: nx2 array of y-coordinate limits for each selected polygon
% =========================================================================
% This function is part of the:
% Geodetic Bayesian Inversion Software (GBIS)
% Software for the Bayesian inversion of geodetic data.
% Copyright: Marco Bagnardi, 2018
%
% Email: gbis.software@gmail.com
%
% Reference: 
% Bagnardi M. & Hooper A, (2018). 
% Inversion of surface deformation data for rapid estimates of source 
% parameters and uncertainties: A Bayesian approach. Geochemistry, 
% Geophysics, Geosystems, 19. https://doi.org/10.1029/2018GC007585
%
% The function may include third party software.
% =========================================================================
% Last update: 8 August, 2018


%% If less than 3 points set xlim/ylim to 0
if(size(xy,1)<3)
    xlims=0;
    ylims=0;
end

%% Initialise variables;
if (~exist('startlevel','var'))
    startlevel = 1;
end
nb        = 0;
err       = 0;
npoints   = 0;
centers   = 0;
polys     = 0;
values    = 0;
%% Open new figure window
if startlevel == 1
    figure('Position', [1, 1, 700, 700]);
end

% If less than 3 points stop quadtree
if(size(xy,1)<3)
    return
end

% Check if patch level has reached maximum level
if(startlevel > maxlevel )
    % Find the convex hull for current set of points
    K = convhull(xy(:,2), xy(:,3));
    plot(xy(K,2), xy(K,3), 'color',[0.8 0.8 0.8]); colorbar; hold on
    return
end

% Determine patch size, centre, and limits
dx = (max(xy(:,2))-min(xy(:,2)));
dy = (max(xy(:,3))-min(xy(:,3)));
cx = min(xy(:,2))+dx/2;
cy = min(xy(:,3))+dy/2;
lim    = [min(xy(:,2)) max(xy(:,2)) min(xy(:,3)) max(xy(:,3))];
xlims  = [lim(1) lim(1) lim(2) lim(2) lim(1)];
ylims  = [lim(3) lim(4) lim(4) lim(3) lim(3)];

% Calculate mean value and variance within cell
mpoly  = mean(val);
spoly  = var(val);

%% If variance < threshold stop quadtree
if( size(xy,1)<100 || spoly < thresh )
    % Find the convex hull for current set of points
    if startlevel == 1
        error(['Quadtree threshold variance (',num2str(thresh),') is larger than data variance (',num2str(spoly),'). Please, lower it.'])
    end
%     if (xy(3,3)==xy(1,3))
%         xy(3,3)=xy(3,3)+0.000002;
%     end
%     if (xy(3,2)==xy(1,2))
%         xy(3,2)=xy(3,2)+0.000002;
%         
%     end
%     save quad_test
    K = convhull(xy(:,2), xy(:,3));
    patch(xy(K,2),xy(K,3),mpoly); hold on;colorbar;
    nb = 1;
    err = sum(abs(val-mpoly));
    npoints = length(val);
    centers = [cx cy];
    values  = mpoly;
    polys= cell(1);
    polys{1}= xy(K,:);
    return
end

%% If variance > threshold continue quadtree
first = 1;
for i = 1:4
    % Quadtree subdivision, calculate coordinates of vertex
    switch(i)
        case 1
            xyv =[cx-dx, cx-dx,  cx,    cx; ...
                cy-dy, cy,     cy,    cy-dy]';
        case 2
            xyv =[cx,    cx,     cx+dx, cx+dx; ...
                cy-dy, cy,     cy,    cy-dy]';
        case 3
            xyv =[cx-dx, cx-dx,  cx,    cx; ...
                cy,    cy+dy,  cy+dy, cy]';
        case 4
            xyv =[cx,    cx,     cx+dx,  cx+dx; ...
                cy,    cy+dy,  cy+dy,  cy]';
    end
    
    % Find ID of points within cell
    if(size(xyv)>0)
        in = find( xy(:,2) <= xyv(3,1) & xy(:,2) >= xyv(1,1) & xy(:,3) <= xyv(2,2) & xy(:,3) >= xyv(1,2) );
        % If there are points continue quadtree (calling same function!).
        if(size(in)>0)
            [pnb,perr,pnpoints,pcenter,pvalue,ppoly,pxlims,pylims] = quadtree(xy(in,:), val(in), thresh, maxlevel, startlevel+1);
            if (pnb>0)
                % Store the results
                nb = nb + pnb;
                err = err + perr;
                if (first == 1)
                    centers = pcenter;
                    values  = pvalue;
                    polys   = ppoly;
                    npoints = pnpoints;
                    xlims   = pxlims;
                    ylims   = pylims;
                    first   = 2;
                else
                    centers = [centers; pcenter];
                    values  = [values, pvalue];
                    polys   = [polys; ppoly];
                    npoints = [npoints; pnpoints];
                    xlims   = [xlims; pxlims];
                    ylims   = [ylims; pylims];
                end
            end
        end
    end
end
return

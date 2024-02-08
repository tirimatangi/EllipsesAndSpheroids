function output = plotEllipses(sharedFocalPoint, otherFocalPoints, diameters, intersections, titleText)
%plotEllipses Plot ellipses with the given focal points and diameters. Also show the points of intersections.
% Inputs:
%   sharedFocalPoint - Coordinate of shared focal points of every ellipse (2x1)
%   otherFocalPoints - Coordinates of N other focal points (2xN)
%   diameters - Diameters of the ellipses along the major axis (Nx1)
%   intersections - Coordinates of M intersections (2xM)

    if nargin < 5
        titleText = "";
    end

    numEllipses = size(otherFocalPoints, 2);
    numPoints = 100;
    ellipseCoord = zeros(2,numPoints, numEllipses);
    theta = linspace(0, 2*pi, numPoints);
    diffCoord = otherFocalPoints - sharedFocalPoint;
    for ii = 1:numEllipses
        a = diameters(ii)/2;  % Major axis
        b = sqrt(a*a - norm(diffCoord(:,ii)/2)^2); % Minor axis
        phase = diffCoord(1,ii) + 1i*diffCoord(2,ii);
        phase = phase / norm(phase); % Rotation of the ellipse
        center = (otherFocalPoints(:,ii) + sharedFocalPoint) / 2; % Center of the ellipse
        for jj = 1:numPoints
            coord = (a * cos(theta(jj))) + 1i*(b * sin(theta(jj))); % axis aligned
            coord = coord * phase;  % rotate
            coord = coord + (center(1) + 1i*center(2));
            ellipseCoord(:,jj,ii) = [real(coord); imag(coord)];
        end
    end
    rgb=['g','r','b'];
    limX = [floor(min(ellipseCoord(1,:))) ceil(max(ellipseCoord(1,:)))];
    limY = [floor(min(ellipseCoord(2,:))) ceil(max(ellipseCoord(2,:)))];
    legendArray = {};

    figure;
    hold on
    % Plot the shared focal point
    plot(sharedFocalPoint(1),sharedFocalPoint(2),'Marker','o','MarkerSize',10,'Color','k');
    legendArray{end+1} = "Shared focal point";

    % Plot non-shared focal points
    for ii=1:numEllipses
        plot(otherFocalPoints(1,ii),otherFocalPoints(2,ii),'Marker','*','MarkerSize',10,'Color',rgb(mod(ii,3)+1));
        legendArray{end+1} = ["focal point #"  num2str(ii)];
    end

    % Plot ellipses
    for ii=1:numEllipses
        plot(ellipseCoord(1,:,ii), ellipseCoord(2,:,ii),'-','Color',rgb(mod(ii,3)+1));
        legendArray{end+1} = ["Ellipse #"  num2str(ii)];
    end


    % Plot the intersections
    for jj=1:size(intersections,2)
        if mod(jj,2) == 0
            mkr = 'pentagram';
        else
            mkr = 'hexagram';
        end
        plot(intersections(1,jj), intersections(2,jj),'Marker',mkr,'MarkerSize',12,'Color','k','LineStyle',':');
        legendArray{end+1} = ["Intersection #" num2str(jj)];
    end

    for jj=1:size(intersections,2)
        for ii=1:numEllipses
            plot([intersections(1,jj) otherFocalPoints(1,ii)], [intersections(2,jj) otherFocalPoints(2,ii)],':','Color',rgb(mod(ii,3)+1))
            plot([intersections(1,jj) sharedFocalPoint(1)], [intersections(2,jj) sharedFocalPoint(2)],':','Color',rgb(mod(ii,3)+1))
        end
    end

    grid on
    axis equal
    ax = gca;
    set (ax, "XTick", limX(1):limX(2));
    set (ax, "YTick", limY(1):limY(2));
    title(titleText)
    legend(legendArray);

    output = gca;
end

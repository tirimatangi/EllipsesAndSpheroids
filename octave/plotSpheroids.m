function output = plotSpheroids(sharedFocalPoint, otherFocalPoints, diameters, intersections, titleText)
%plotEllipses Plot spheroids with the given focal points and diameters. Also show the points of intersections.
% Inputs:
%   sharedFocalPoint - Coordinate of shared focal points of every ellipse (2x1)
%   otherFocalPoints - Coordinates of N other focal points (2xN)
%   diameters - Diameters of the ellipses along the major axis (Nx1)
%   intersections - Coordinates of M intersections (2xM)

    if nargin < 5
        titleText = "";
    end

    numPoints = 50;
    numSpheroids = size(otherFocalPoints, 2);
    spheroidCoord = zeros(3,numPoints, numSpheroids);
    theta = linspace(0, 2*pi, numPoints);
    diffCoord = otherFocalPoints - sharedFocalPoint;

    figure;
    hold on
    grid on
    axis equal
    legendArray = {};

    for ii = 1:numSpheroids
        a = diameters(ii)/2;  % Major axis
        b = sqrt(a*a - norm(diffCoord(:,ii)/2)^2); % Minor axis
        helperEllipse = zeros(3, numPoints*numPoints);
        helperIndex = 0;
        for jj = 1:numPoints
            coord = [a * cos(theta(jj)) , b * sin(theta(jj))];
            for kk = 1:numPoints
                helperIndex = helperIndex + 1;
                helperEllipse(:,helperIndex) = [coord(1); cos(theta(kk))*coord(2); sin(theta(kk))*coord(2)];
            end
        end

        center = (otherFocalPoints(:,ii) + sharedFocalPoint) / 2; % Center of the spheroid

        % Cartesian x axis will be mapped into this vector,
        % which is the direction of the major axis.
        Vaxis = diffCoord(:,ii);
        Vaxis = Vaxis / norm(Vaxis);
        % Decide another axis orthogonal the Xaxis
        [~,minIndex] = min(abs(axis));
        switch minIndex
            case 1
                Uaxis = [1;0;0];
            case 2
                Uaxis = [0;1;0];
            case 3
                Uaxis = [0;0;1];
        end
        Waxis = Uaxis - (Uaxis' * Vaxis) * Vaxis;
        Waxis = Waxis / norm(Waxis);
        Uaxis=cross(Waxis, Vaxis);
        % This matrix rotates the helper ellipse into place
        Rot = [Vaxis Waxis Uaxis];

        spheroid = Rot * helperEllipse + center;

        plot3(spheroid(1,:), spheroid(2,:), spheroid(3,:),'.','MarkerSize',10);
        legendArray{end+1} = ["Spheroid " num2str(ii)];
        % Plot Rx->Tx pairs
        plot3([sharedFocalPoint(1), otherFocalPoints(1,ii)], ...
            [sharedFocalPoint(2), otherFocalPoints(2,ii)], ...
            [sharedFocalPoint(3), otherFocalPoints(3,ii)],'x-','LineWidth',3);
        legendArray{end+1} = ["Foci axis " num2str(ii)];
    end

    % Plot the intersections
    numIntersections = size(intersections,2);
    if numIntersections == 1
            plot3(intersections(1), intersections(2), intersections(3), ...
                'ro', 'MarkerSize',16,'MarkerFaceColor','auto','LineWidth',2);
            legendArray{end+1} = "Intersection";
    elseif numIntersections > 1
        for jj=1:size(intersections,2)
            if mod(jj,2) == 0
                mkr = 'pentagram';
            else
                mkr = 'hexagram';
            end
            plot3(intersections(1,jj), intersections(2,jj), intersections(3,jj), ...
                'Color','r','Marker',mkr, 'MarkerSize',16,'MarkerFaceColor','auto','LineWidth',2);
            legendArray{end+1} = ["Intersection #"  num2str(jj)];
        end
    end

    xlabel("XX")
    ylabel("YY")
    zlabel("ZZ")
    legend(legendArray)
    title(titleText)
    output = gca;
end

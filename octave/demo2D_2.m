%% This test demonstrates intersections of 2 ellipses

% Coordinate of the shared focal point
sharedFocalPoint = [-10; -20];

% Coordinates of the other focal points of 2 ellipses.
numEllipses = 2;
otherFocalPoints = zeros(2, numEllipses);

%% (Case 1: The major axes are collinear and upright)
otherFocalPoints(:,1) = [0; 20] + sharedFocalPoint;
otherFocalPoints(:,2) = [0; 10] + sharedFocalPoint;

% Select point p at which the ellipses intersect and determine the diameters accordingly.
% There will be 2 points of intersection, one of which will be p.
p = [4; 3] + sharedFocalPoint;
diameters = threePointDistance(otherFocalPoints, p, sharedFocalPoint)';

% Solve intersections and an estimate on the model error.
% The bigger the error, the further away the estimated intersections are from the real ones. 
% The model error is not exacly the actual distance between 
% the estimate and the correct intersection so it should be taken with a pinch of salt.
[intersections, relativeModelError] = solveEllipseIntersections(sharedFocalPoint, otherFocalPoints, diameters);

assert(min(vecnorm(intersections - p)) < 1e-6, "Intersection point is not correct")
assert(max(relativeModelError) < 1e-6, "Model error in the solver is too large")

% Calculate ellipse coordinates and plot the coordinates and the intersection points.
plotEllipses(sharedFocalPoint, otherFocalPoints, diameters, intersections, "2 ellipses, Case 1");

disp("Case 1: intersection coordinates as columns: ");
disp(intersections);

%% (Case 2: The major axes are collinear and tilted)
otherFocalPoints(:,1) = [-10; 20] + sharedFocalPoint;
otherFocalPoints(:,2) = [-5; 10] + sharedFocalPoint;

% Select point p at which the ellipses intersect and determine the diameters accordingly.
% There will be 2 points of intersection, one of which will be p.
p = [4; 3] + sharedFocalPoint;
diameters = threePointDistance(otherFocalPoints, p, sharedFocalPoint)';

% Solve intersections and an estimate on the model error.
% The bigger the error, the further away the estimated intersections are from the real ones. 
% The model error is not exacly the actual distance between 
% the estimate and the correct intersection so it should be taken with a pinch of salt.
[intersections, relativeModelError] = solveEllipseIntersections(sharedFocalPoint, otherFocalPoints, diameters);

assert(min(vecnorm(intersections - p)) < 1e-6, "Intersection point is not correct")
assert(max(relativeModelError) < 1e-6, "Model error in the solver is too large")

% Calculate ellipse coordinates and plot the coordinates and the intersection points.
plotEllipses(sharedFocalPoint, otherFocalPoints, diameters, intersections, "2 ellipses, Case 2");

disp("Case 2: intersection coordinates as columns:");
disp(intersections);

%% (Case 3: General case with 2 ellipses)
otherFocalPoints(:,1) = [10; 4] + sharedFocalPoint;
otherFocalPoints(:,2) = [0; 10] + sharedFocalPoint;

% Select point p at which the ellipses intersect and determine the diameters accordingly.
% There will be 2 points of intersection, one of which will be p.
p = [4; 3] + sharedFocalPoint;
diameters = threePointDistance(otherFocalPoints, p, sharedFocalPoint)';

% Solve intersections and an estimate on the model error.
% The bigger the error, the further away the estimated intersections are from the real ones. 
% The model error is not exacly the actual distance between 
% the estimate and the correct intersection so it should be taken with a pinch of salt.
[intersections, relativeModelError] = solveEllipseIntersections(sharedFocalPoint, otherFocalPoints, diameters);

assert(min(vecnorm(intersections - p)) < 1e-6, "Intersection point is not correct")
assert(max(relativeModelError) < 1e-6, "Model error in the solver is too large")

% Calculate ellipse coordinates and plot the coordinates and the intersection points.
plotEllipses(sharedFocalPoint, otherFocalPoints, diameters, intersections, "2 ellipses, Case 3");

disp("Case 3: intersection coordinates as columns:");
disp(intersections);
%% This test demonstrates intersections of 3 ellipses

% Coordinate of the shared focal point
sharedFocalPoint = [-10; -20];

% Coordinates of the other focal points of 2 ellipses.
numEllipses = 3;
otherFocalPoints = zeros(2, numEllipses);

%% (Case 1: The major axes are collinear and upright)
otherFocalPoints(:,1) = [0; 20] + sharedFocalPoint;
otherFocalPoints(:,2) = [0; 10] + sharedFocalPoint;
otherFocalPoints(:,3) = [0; -16] + sharedFocalPoint;

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
plotEllipses(sharedFocalPoint, otherFocalPoints, diameters, intersections, "3 ellipses, Case 1");

disp("Case 1: intersection coordinates as columns: ");
disp(intersections);

%% (Case 2: The major axes are collinear and tilted)
otherFocalPoints(:,1) = [-10; 20] + sharedFocalPoint;
otherFocalPoints(:,2) = [-5; 10] + sharedFocalPoint;
otherFocalPoints(:,3) = [-8; 16] + sharedFocalPoint;

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
plotEllipses(sharedFocalPoint, otherFocalPoints, diameters, intersections, "3 ellipses, Case 2");

disp("Case 2: intersection coordinates as columns:");
disp(intersections);

%% (Case 3: General case with 3 ellipses)
otherFocalPoints(:,1) = [10; 4] + sharedFocalPoint;
otherFocalPoints(:,2) = [0; 10] + sharedFocalPoint;
otherFocalPoints(:,3) = [10; 0] + sharedFocalPoint;

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
plotEllipses(sharedFocalPoint, otherFocalPoints, diameters, intersections, "3 ellipses, Case 3");

disp("Case 3: (One intersection only) Intersection coordinates as columns:");
disp(intersections);

%% (Case 4: Non-solvable case with 3 ellipses)
% This is an example of the solver giving a wrong answer and returning a high
% model error when the requirement that the ellipses intersect at one point does not hold.
otherFocalPoints(:,1) = [10; 4] + sharedFocalPoint;
otherFocalPoints(:,2) = [0; 10] + sharedFocalPoint;
otherFocalPoints(:,3) = [10; 0] + sharedFocalPoint;

% Select point p at which the ellipses intersect and determine the diameters accordingly.
% There will be 2 points of intersection, one of which will be p.
p = [4; 3] + sharedFocalPoint;
diameters = threePointDistance(otherFocalPoints, p, sharedFocalPoint)';

% The 3rd ellipse does not intersect at the same point where the 1st and 2nd ellipses intersect
pp = p + [1; 1];
diameters(1) = threePointDistance(otherFocalPoints(:,1), pp, sharedFocalPoint)';

% Solve intersections and an estimate on the model error.
% The bigger the error, the further away the estimated intersections are from the real ones.
% The model error is not exacly the actual distance between
% the estimate and the correct intersection so it should be taken with a pinch of salt.
[intersections, relativeModelError] = solveEllipseIntersections(sharedFocalPoint, otherFocalPoints, diameters);

% Calculate ellipse coordinates and plot the coordinates and the intersection points.
plotEllipses(sharedFocalPoint, otherFocalPoints, diameters, intersections, ...
             ["3 ellipses intersecting at different points, rel.model error=" num2str(relativeModelError)]);

disp(["Case 4: intersection with model error " num2str(relativeModelError)]);
disp(intersections);

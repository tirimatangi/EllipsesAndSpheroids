%% This test demonstrates intersections of 3 spheroids

% Coordinate of the shared focal point
sharedFocalPoint = [-10; -20; 30];

% Coordinates of the other focal points of 2 spheroids.
numSpheroids = 3;
otherFocalPoints = zeros(3, numSpheroids);

%% (Case 1: The major axes are on the same yz-plane where x is constant)
otherFocalPoints(:,1) = [0; 20; 1] + sharedFocalPoint;
otherFocalPoints(:,2) = [0; 10; -2] + sharedFocalPoint;
otherFocalPoints(:,3) = [0; -16; 3] + sharedFocalPoint;

% Select point p at which the spheroids intersect and determine the diameters accordingly.
% There will be 2 points of intersection, one of which will be p.
p = [4; 3; 11] + sharedFocalPoint;
diameters = threePointDistance(otherFocalPoints, p, sharedFocalPoint)';

% Solve intersections and an estimate on the model error.
% The bigger the error, the further away the estimated intersections are from the real ones. 
% The model error is not exacly the actual distance between 
% the estimate and the correct intersection so it should be taken with a pinch of salt.
[intersections, relativeModelError] = solveEllipseIntersections(sharedFocalPoint, otherFocalPoints, diameters);

assert(min(vecnorm(intersections - p)) < 1e-6, "Intersection point is not correct")
assert(max(relativeModelError) < 1e-6, "Model error in the solver is too large")

% Calculate spheroid coordinates and plot the coordinates and the intersection points.
plotSpheroids(sharedFocalPoint, otherFocalPoints, diameters, intersections, "3 spheroids, Case 1");

disp("Case 1: intersection coordinates as columns: ");
disp(intersections);

%% (Case 2: The major axes are on the same tilted plane (x and y coordinates are linearly dependent))
otherFocalPoints(:,1) = [10; 20; 1] + sharedFocalPoint;
otherFocalPoints(:,2) = [5; 10; -2] + sharedFocalPoint;
otherFocalPoints(:,3) = [-8; -16; 3] + sharedFocalPoint;

% Select point p at which the spheroids intersect and determine the diameters accordingly.
% There will be 2 points of intersection, one of which will be p.
p = [4; 3; 11] + sharedFocalPoint;
diameters = threePointDistance(otherFocalPoints, p, sharedFocalPoint)';

% Solve intersections and an estimate on the model error.
% The bigger the error, the further away the estimated intersections are from the real ones. 
% The model error is not exacly the actual distance between 
% the estimate and the correct intersection so it should be taken with a pinch of salt.
[intersections, relativeModelError] = solveEllipseIntersections(sharedFocalPoint, otherFocalPoints, diameters);

assert(min(vecnorm(intersections - p)) < 1e-6, "Intersection point is not correct")
assert(max(relativeModelError) < 1e-6, "Model error in the solver is too large")

% Calculate spheroid coordinates and plot the coordinates and the intersection points.
plotSpheroids(sharedFocalPoint, otherFocalPoints, diameters, intersections, "3 spheroids, Case 2");

disp("Case 2: intersection coordinates as columns:");
disp(intersections);

%% (Case 3: z coordinates of the major axes are linearly dependent on x- and y-coordinates.)
otherFocalPoints(:,1) =[10; 1; 10+1] + sharedFocalPoint;
otherFocalPoints(:,2) = [5; -2; 5-2] + sharedFocalPoint;
otherFocalPoints(:,3) = [-8; 3; -8+3] + sharedFocalPoint;

% Select point p at which the spheroids intersect and determine the diameters accordingly.
% There will be 2 points of intersection, one of which will be p.
p = [4; 3; 11] + sharedFocalPoint;
diameters = threePointDistance(otherFocalPoints, p, sharedFocalPoint)';

% Solve intersections and an estimate on the model error.
% The bigger the error, the further away the estimated intersections are from the real ones. 
% The model error is not exacly the actual distance between 
% the estimate and the correct intersection so it should be taken with a pinch of salt.
[intersections, relativeModelError] = solveEllipseIntersections(sharedFocalPoint, otherFocalPoints, diameters);

assert(min(vecnorm(intersections - p)) < 1e-6, "Intersection point is not correct")
assert(max(relativeModelError) < 1e-6, "Model error in the solver is too large")

% Calculate spheroid coordinates and plot the coordinates and the intersection points.
plotSpheroids(sharedFocalPoint, otherFocalPoints, diameters, intersections, "3 spheroids, Case 3");

disp("Case 3: intersection coordinates as columns:");
disp(intersections);
%% (Case 4: General case with 3 spheroids)
otherFocalPoints(:,1) = [-2; -5; 1] + sharedFocalPoint;
otherFocalPoints(:,2) = [0; 10; 2] + sharedFocalPoint;
otherFocalPoints(:,3) = [10; 0; 3] + sharedFocalPoint;

% Select point p at which the spheroids intersect and determine the diameters accordingly.
% There will be 2 points of intersection, one of which will be p.
p = [4; 3; 11] + sharedFocalPoint;
diameters = threePointDistance(otherFocalPoints, p, sharedFocalPoint)';

% Solve intersections and an estimate on the model error.
% The bigger the error, the further away the estimated intersections are from the real ones. 
% The model error is not exacly the actual distance between 
% the estimate and the correct intersection so it should be taken with a pinch of salt.
[intersections, relativeModelError] = solveEllipseIntersections(sharedFocalPoint, otherFocalPoints, diameters);

assert(min(vecnorm(intersections - p)) < 1e-6, "Intersection point is not correct")
assert(max(relativeModelError) < 1e-6, "Model error in the solver is too large")

% Calculate spheroid coordinates and plot the coordinates and the intersection points.
plotSpheroids(sharedFocalPoint, otherFocalPoints, diameters, intersections, "3 spheroids, Case 4");

disp("Case 4: intersection coordinates as columns:");
disp(intersections);
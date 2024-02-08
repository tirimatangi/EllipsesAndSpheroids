function distance = threePointDistance(point1Array, point2, point3)
%measuredDistances Calculate distance from each Point1 to Point2 and further on to Point3 
% Inputs:
%   point1Array - array of Point1 coordinates, DIM x N array
%   point2 -  DIM x 1 array
%   point3 - DIM x 1 array
% Output:
%  distances - point1(k) -> point2 -> point3 distance, 1 x N array
    distance = vecnorm(point1Array - point2, 2, 1);
    distance = distance + norm(point2 - point3);
end
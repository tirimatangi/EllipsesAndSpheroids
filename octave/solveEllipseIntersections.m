function [intersectionPoints, relativeModelError] = solveEllipseIntersections(sharedFocalPoint, otherFocalPoints, diameters)
%solveIntersections Calculates intersections of ellipses and spheroids with a shared focal point.
%   Calculates intersections of ellipses (2D) or spheroids (3D). All ellipses are assumed to have the first focal point  
%   at the same position and the second focal point at the given coordinates.
%   The algorithm for 3D spheroids is described below. The 2D case is similar but simpler.
%
%   The diameters along the major semiaxes determine the sum of distances from the first focal point to 
%   any point on the locus of the ellipse or spheroid and further on to the second focal point.
%   Without loss of generality we assume that the first focal point of each spheroid is at the origin 
%   and fp = otherFocalPoints - sharedFocalPoint is the set of the second focal points.
%   Now assume that [x y z] is an intersection of the spheroids. Its distance from the origin is denoted by w.
%   Now the intersection point must satisfy the following N+1 conditions where N is the number of intersecting spheroids:
%
%                                           x^2 + y^2 + z^2 = w^2
%       (x - fp(1,1))^2 + (y - fp(2,1))^2 + (z - fp(3,1))^2 = (diam(1) - w)^2
%       ...
%       (x - fp(1,N))^2 + (y - fp(2,N))^2 + (z - fp(3,N))^2 = (diam(N) - w)^2
%
% By expanding the squares in N last equations and subtracting the first equation
% from the expanded squares, we get an N-by-4 linear system  T*[x y z w]' = r as follows:
%
%   [-2 fp(1,1) -2 fp(2,1) -2 fp(3,1) diam(1)] [x] = [diam(1)^2 - fp(1,1)^2 - fp(2,1)^2 - fp(3,1)^2]
%   [                   ...                  ] [y] = [      ....                       ]
%   [                   ...                  ] [z] = [      ....                       ]
%   [-2 fp(1,N) -2 fp(2,N) -2 fp(3,N) diam(N)] [w] = [diam(N)^2 - fp(1,N)^2 - fp(2,N)^2 - fp(3,N)^2]
%
%   So T = -2 * [fp(1,:)' fp(2,:)' fp(3,:)' diam] and r =  diam.^2 - fp(1,:)'.^2 - fp(2,:)'.^2 - fp(3,:)'.^2;
%
%   Notice that a row of T may be linearly dependent on the other rows. Those cases can be dealt with
%   by writing the dependent row as a linear combination of the other rows and using the fact that
%   x^2 + y^2 + z^2 = w^2 to solve the coefficients of the linear combination. This will be explained in
%   comments in the code below.
%
%   Inputs:
%       sharedFocalPoint - The focal point which is shared by all ellipses / spheroids, DIM-by-1 vector
%       otherFocalPoints - The other (non-shared) focal points as DIM-by-N matrix.
%          In 2D case, N should be 2 or 3. In 3D case, N should be 3 or 4.
%       diameters - Diameters of the ellipses along the major semiaxes as an N-by-1 vector.
%   Output:
%       intersectionPoints - Coordinates of the points where all N ellipses intersect. 
%                            A 2-by-M or 3-by-M matrix where M is 1 of the intersection point is unique and 2 if there are two candidates.
%                            If the matrix is empty, there is no common intersection point.
%                            Note that if there is noise or other inaccuracy in the set of diameters, the all ellipses/spheroids
%                            may not actually intersect at the given points. In this case distanceError will be large
%                            so the inaccuracy can be detected.
%       relativeModelError - Relative difference between the actual squared distance of the calculated intersection point(s)
%                            from the shared focal point (i.e. x^2+y^2+z^2) and the squared distance predicted by
%                            the model (i.e. w^2). So relativeModelError = sqrt(abs(x^2+y^2+z^2 - w^2) / x^2+y^2+z^2).
%                            The larger the value, the less reliable the result is. 
%                            The value is very near zero if all ellipses/spheroids intersect at the same points.

assert(nargin == 3)
% 2D or 3D mode?

dimension = size(otherFocalPoints, 1);
numFoci = size(otherFocalPoints, 2);

assert(numFoci == length(diameters), "The number of focal points and the number of diameters must be the same.");
assert(size(diameters, 2) == 1, "Diameters must be an N-by-1 vector")
assert(size(sharedFocalPoint,1) == dimension)
assert(size(sharedFocalPoint,2) == 1)

% Normalize the coordinates so that the shared focal point is at the origin.
focalPoints = otherFocalPoints - sharedFocalPoint;
switch dimension
    case 2
        assert(numFoci == 2 || numFoci == 3, "In 2D case the number of points should be 2 or 3")
         intersectionPoints = solveEllipseIntersections2D(focalPoints, diameters);
    case 3
        assert(numFoci == 3 || numFoci == 4, "In 3D case the number of points should be 3 or 4")
        intersectionPoints = solveSpheroidIntersections3D(focalPoints, diameters);
    otherwise
        error("The dimension of focal points must be 2 or 3")
end


if isempty(intersectionPoints)
    relativeModelError = Inf;
    return
end

% Sometimes the imaginary components can be non-zero due to numerical errors in the solver.
if ~isreal(intersectionPoints)
    % The points are conjugates, so only one real part is enough.
    intersectionPoints = real(intersectionPoints(:,1));
end
if nargout > 1
    w = intersectionPoints(dimension+1,:);
    xyz = intersectionPoints(1:dimension,:);
    relativeModelError =  sqrt(abs(vecnorm(xyz).^2 - w.^2)) ./ (vecnorm(xyz) + realmin);
end

% Undo the normalization
intersectionPoints = intersectionPoints(1:dimension,:) + sharedFocalPoint;
end


function intersectionPoints = solveEllipseIntersections2D(fp, diam)
%   Inputs:
%       fp - Second focal points as a 2-by-N matrix (the first focal points are always at the origin) 
%       diam - Diameters of the ellipses along the major semiaxes, N-by-1 vector.
%   Output:
%       intersectionPoints - Coordinates of the points where all N ellipses intersect, or an empty vector if they don't intersect.
%
    DIM = 2;
    intersectionPoints = [];
    numFoci = length(diam); % Number of ellipses (must be DIM or DIM+1)

%  An intersection point (x,y) of the spheroids and it's distance w = sqrt(x^2 + y^2) from the origin 
% (i.e. the second focal point) satifies the following conditions:
%       (x - fp(1,1))^2 + (y - fp(2,1))^2 = (diam(1) - w)^2
%       ...
%       (x - fp(1,N))^2 + (y - fp(2,N))^2 = (diam(N) - w)^2
%                               x^2 + y^2 = w^2
% By expanding the squares in N first equations and using the last equation, 
% we get a linear problem T * [x y w]' = r  with matrix T and vector r defined as as follows:
    T = -2 * [fp(1,:)' fp(2,:)' diam];
    r = diam.^2 - fp(1,:)'.^2 - fp(2,:)'.^2;
    
    % QR transform will later help with sorting out corner cases later
    [~,R] = qr(T);
    
    % Are there tiny elements on the diagonal of R? 
    % If so, it means that some columns of A are linearly dependent.
    zeroDiagIndex = findSmallIndices(diag(R));

    % If there are fewer than DIM+1 ellipses, "add a fictional zero row to the end of R",´.
    % This means that the (fictional) last element of the diagonal is zero 
    % (i.e. it's index belongs in the list of small elements).
    % This will trigger the correct solution option later. Don't worry, no panic, read on.
    if isempty(zeroDiagIndex) &&  numFoci <= DIM
        zeroDiagIndex = DIM + 1;
    end
    
    % If there are more than one zero on the diagonal, the intersections can not be solved.
    if length(zeroDiagIndex) > 1
        return 
    end

    % If there are no zeros on the diagonal, the problem can be solved
    if isempty(zeroDiagIndex)
        xyw = T \ r; % == R \ (Q' * r);
        intersectionPoints = xyw;
        return
    end

    % Now matrix T is made of 3 columns: T = [A B C] and the problem to be solved is
    % T * [x y w]' = Ax + By + Cw = r  such that  x^2 + y^2 = w^2
    % Now one of the diagonals of the upper diagonal matrix R is zero,
    % meaning that the corresponding column of T can be expressed as a linear combination
    % of the preceeding columns. 
    % Let's go over all possible positions of the zero diagonal and see what it means.
    switch zeroDiagIndex
        case 1
            % R(1,1) = 0, meaning that column A = 0. 
            % Solve By + Cw = r and use constraint x^2 = w^2 - y^2
            yw = T(:,2:3)\r; % Now yw(1) = y, yw(2) = w
            x = sqrt(yw(2)^2 - yw(1)^2);
            x = [x  -x];
            y = [yw(1) yw(1)];
            w = [yw(2) yw(2)]; % Uncomment if w is needed
        case 2
            % R(2,2) = 0, meaning that column B = g * A for some g.
            % Let's find what g is:
            g = R(1,2) / R(1,1);
            % Now we have  Ax + (g * A)y + Cw = r  such that  x^2 + y^2 = w^2
            % --> A (x + g*y) + Cw = r. 
            % Denote h=x + g*y, meaning x = h - g*y.
            % We can now solve h and w without using column B:
            hw = T(:,[1 3]) \ r; % Now hw(1) = h, hw(2) = w
            % Now w^2 = x^2 + y^2 = (h - g*y)^2 + y^2
            % --> (g^2 + 1)*y^2 - 2*g*h*y + (h^2 - w^2) = 0.
            %  This is a 2nd degree polynomial on y. Let's solve it.
            y = poly2zeros(g^2 + 1, -2*g*hw(1), hw(1)^2 - hw(2)^2);
            % Now that we know y, we also know x = -g*y + h
            x = hw(1) - g*y ;
            w = [hw(2) hw(2)]; % Uncomment if w is needed
        case 3
            % R(3,3) = 0, meaning that column C = g1*A + g2*B for some g1,g2.
            % Let's find what g is:
            g = R(1:2,1:2) \ R(1:2,3);
            % Now A*x + B*y + (g1*A + g2*B)*w = r such that  x^2 + y^2 = w^2
            % --> A*(x+g1*w) + B*(y+g2*w) = r.
            % Denote h1 = x+g1*w, h2 = y+g2*w, meaning x = h1 - g1*w, y = h2 - g2*w
            % Now solve h1 and h2 using columns A and B
            h = T(:,1:2) \ r;
            % Now  w^2 = x^2 + y^2 = (h1 - g1*w)^2 + (h2 - g2*w)^2
            % --> (g1^2 + g2^2 - 1)*w^2 - 2*(g1*h1+g2*h2)*w + (h1^2+h2^2) = 0
            %  This is a 2nd degree polynomial on w. Let's solve it.
            w = poly2zeros(g(1)^2 + g(2)^2 - 1, -2*(g(1)*h(1) + g(2)*h(2)), h(1)^2 + h(2)^2);
            % Now we know w, g and h so we know x and y
            xy = h - g*w;
            x = xy(1,:);
            y = xy(2,:);
        otherwise
            error("small diags should be 1 2 or 3 (" + zeroDiagIndex + ")")
    end
    intersectionPoints = [x;y;w];
end

function intersectionPoints = solveSpheroidIntersections3D(fp, diam)
%   Inputs:
%       fp - Second focal points as a 3-by-N matrix (the first focal points are always at the origin) 
%       diam - Diameters of the ellipses along the major semiaxes, N-by-1 vector.
%   Output:
%       intersectionPoints - Coordinates of the points where all N ellipses intersect, or an empty vector if they don't intersect.
    DIM = 3;
    intersectionPoints = [];
    numFoci = length(diam); % Number of ellipses (must be DIM or DIM+1)

%  An intersection point (x,y,z) of the spheroids and it's distance w = sqrt(x^2 + y^2 + z^2) from the origin 
% (i.e. the second focal point) satifies the following conditions:
%       (x - fp(1,1))^2 + (y - fp(2,1))^2 + (z - fp(3,1))^2 = (diam(1) - w)^2
%       ...
%       (x - fp(1,N))^2 + (y - fp(2,N))^2 + (z - fp(3,N))^2 = (diam(N) - w)^2
%                                           x^2 + y^2 + z^2 = w^2
% By expanding the squares in N first equations and using the last equation, 
% we get a linear problem T * [x y z w]' = r  with matrix T and vector r defined as as follows:
    T = -2 * [fp(1,:)' fp(2,:)' fp(3,:)' diam];
    r = diam.^2 - fp(1,:)'.^2 - fp(2,:)'.^2 - fp(3,:)'.^2;
    
    % QR transform will later help with sorting out corner cases later
    [~,R] = qr(T);
    
    % Are there tiny elements on the diagonal of R? 
    % If so, it means that some columns of A are linearly dependent.
    zeroDiagIndex = findSmallIndices(diag(R));

    % If there are fewer than DIM+1 ellipses, "add a fictional zero row to the end of R",´.
    % This means that the (fictional) last element of the diagonal is zero 
    % (i.e. it's index belongs in the list of small elements).
    % This will trigger the correct solution option later. Don't worry, no panic, read on.
    if isempty(zeroDiagIndex) && numFoci <= DIM
        zeroDiagIndex = DIM + 1;
    end
    
    % If there are more than one zero on the diagonal, the intersections can not be solved.
    if length(zeroDiagIndex) > 1
        return 
    end
    
    % If there are no zeros on the diagonal, the problem can be solved
    if isempty(zeroDiagIndex)
        xyzw = T \ r; % == R \ (Q' * r);
        intersectionPoints = xyzw;
        return
    end
    % Now matrix T is made of 4 columns: T = [A B C D] and the problem to be solved is
    % T * [x y z w]' = Ax + By + Cz + Dw = r  such that  x^2 + y^2 + z^2 = w^2.
    % Now one of the diagonals of the upper diagonal matrix R is zero,
    % meaning that the corresponding column of T can be expressed as a linear combination
    % of the preceeding columns. 
    % Let's go over all possible positions of the zero diagonal and see what it means.
    switch zeroDiagIndex
        case 1
            % R(1,1) = 0, meaning that column A = 0. 
            % Solve By + Cz + Dw = r and use constraint x^2 = w^2 - y^2 - z^2
            yzw = T(:,2:4)\r; % Now yzw(1) = y, yzw(2) = z, yzw(3) = w
            x = sqrt(yzw(3)^2 - yzw(1)^2 - yzw(2)^2);
            x = [x  -x]; % x and -x both are solutions.
            y = [yzw(1) yzw(1)];
            z = [yzw(2) yzw(2)];
            w = [yzw(3) yzw(3)]; % Uncomment if w is needed
        case 2
            % R(2,2) = 0, meaning that column B = g * A for some g.
            % Let's find what g is:
            g = R(1,2) / R(1,1);
            % Now we have  Ax + (g * A)y + Cz + Dw = r  such that  x^2 + y^2 + z^2 = w^2
            % --> A (x + g*y) + Cz + Dw = r. 
            % Denote h=x + g*y, meaning x = -g*y + h.
            % We can now solve h and w without using column B:
            hzw = T(:,[1 3 4]) \ r; % Now hzw(1) = h, hzw(2) = z, hzw(3) = w
            % Now w^2 = x^2 + y^2 + z^2 = (-g*y + h)^2 + y^2 + z^2
            % --> (g^2 + 1)*y^2 - 2*g*h*y + (h^2 + z^2 - w^2) = 0.
            %  This is a 2nd degree polynomial in y. Let's solve it.
            y = poly2zeros(g^2 + 1, -2*g*hzw(1), hzw(1)^2 + hzw(2)^2 - hzw(3)^2);
            % Now that we know y, we also know x = -g*y + h
            x = -g*y + hzw(1);
            z = [hzw(2) hzw(2)];
            w = [hzw(3) hzw(3)]; % Uncomment if w is needed
        case 3
            % R(3,3) = 0, meaning that column C = g1*A + g2*B for some g1,g2.
            % Let's find what g is:
            g = R(1:2,1:2) \ R(1:2,3);
            % Now A*x + B*y + (g1*A + g2*B)*z + Dw = r such that  x^2 + y^2 + z^2 = w^2
            % --> A*(x+g1*z) + B*(y+g2*z) + Dw = r.
            % Denote h1 = x+g1*z, h2 = y+g2*z, meaning x = h1 - g1*z, y = h2 - g2*z
            % Now solve h1 and h2 using columns A, B and D
            hw = T(:,[1 2 4]) \ r; % Now h(1:2) are h1,h2 and hw(3) = w
            % Now w^2 = x^2 + y^2 + z^2 = (h1 - g1*z)^2 + (h2 - g2*z)^2 + z^2
            % --> (g1^2 + g2^2 + 1)*z^2 - 2*(g1*h1 + g2*h2)*z + (h1^2 + h2^2 - w^2) = 0
            % This is a 2nd degree polynomial in z. Let's solve it
            z = poly2zeros(g(1)^2 + g(2)^2 + 1, -2*(g(1)*hw(1) + g(2)*hw(2)), hw(1)^2 + hw(2)^2 - hw(3)^2);
            xy = hw(1:2) - g*z;
            x = [xy(1,1) xy(1,2)];
            y = [xy(2,1) xy(2,2)];
            w = [hw(3) hw(3)];
        case 4
            % R(4,4) = 0, meaning that column D = g1*A + g2*B + g3*C for some g1,g2,g3.
            % Let's find what g is:
            g = R(1:3,1:3) \ R(1:3,4);
            % Now A*x + B*y + C*x + (g1*A + g2*B + g3*C)*w = r such that  x^2 + y^2 +z^2 = w^2
            % --> A*(x+g1*w) + B*(y+g2*w) + C*(z+g3*w) = r.
            % Denote h1 = x+g1*w, h2 = y+g2*w, h3 = z+g3*w, 
            % meaning x = h1 - g1*w, y = h2 - g2*w, z = h3 - g3*w
            % Now solve h1, h2, h3 using columns A, B and C
            h = T(:,1:3) \ r;

            % Now w^2 = x^2 + y^2 + z^2 = (h1 - g1*w)^2 + (h2 - g2*w)^2 + (h3 - g3*w)^2
            % --> (g1^2 + g2^2 + g3^2 - 1)*w^2 - 2*(g1*h1+g2*h2+g3*h3)*w + (h1^2+h2^2+h3^2) = 0
            %  This is a 2nd degree polynomial on w. Let's solve it.
            w = poly2zeros(g(1)^2 + g(2)^2 + g(3)^2 - 1, ...
                           -2*(g(1)*h(1)+g(2)*h(2)+g(3)*h(3)), ...
                           h(1)^2 + h(2)^2 + h(3)^2);
            % Now we know w, g and h so we know x, y and z
            xyz = h - g*w;
            x = xyz(1,:);
            y = xyz(2,:);
            z = xyz(3,:);
        otherwise
            error("small diags should be 1 2 or 3 (" + zeroDiagIndex + ")")
    end
    intersectionPoints = [x;y;z;w];
end

function z = poly2zeros(a,b,c)
% poly2zeros  Finds zeros of polynomial a*x^2 + b*x + c
    determinant = sqrt(b*b - 4*a*c);
    z = (-b + [determinant -determinant]) / (2*a);
end

function smallIndices = findSmallIndices(v)
% smallIndidices Finds "very small" indices of vector v
  smallIndices = find(abs(v) <= sqrt(eps(norm(v))));
end

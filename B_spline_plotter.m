clear all
clc

% Main script for generating and plotting a B-spline curve

%% Control Point Coordinates
controlPointsX = [0 2 4 5 7 9 10];
controlPointsY = [1 6 8 2 8 7 4];

degree = 3; % Degree of the B-spline (degree = k - 1)
k = degree + 1;
numControlPoints = length(controlPointsX) - 1; % Number of control points is n+1

%% Function Call: Compute and plot B-spline curve
bsplineCurve(numControlPoints, k, controlPointsX, controlPointsY)

%% Plot control points
stem(controlPointsX, controlPointsY, 'k', 'LineStyle', 'none', 'Marker', '.', ...
    'MarkerSize', 12, 'DisplayName', 'Control Points')
title('B-spline Curve')
xlabel('X')
ylabel('Y')
legend
grid on
axis equal
hold off

%% Functions

% Function to compute and plot a B-spline curve
function bsplineCurve(n, k, x, y)
    % Inputs:
    % n - number of control points
    % k - order of the B-spline
    % x, y - control point coordinates

    % Generate the knot vector
    knotVector = generateKnotVector(n, k);

    % Define parametric range for t
    tMin = min(knotVector);
    tMax = max(knotVector);
    numPoints = 1000; % Number of points to compute on the curve
    t = linspace(tMin, tMax, numPoints);

    % Initialize storage for the curve points
    curvePoints = zeros(2, length(t));

    % Compute the B-spline curve
    for j = 1:length(t)
        px = 0; py = 0;
        
        for i = 1:(n+1)
            basisValue = basisFunction(i, k, t(j), knotVector);  % Basis function N(i,k) at t(j)
            px = px + x(i) * basisValue;
            py = py + y(i) * basisValue;
        end
        
        curvePoints(:, j) = [px; py];  % Store the computed curve points
    end
    
    % Final Curve points
    curvePoints(:, end) = [x(end); y(end)]; % This is to define the curve at tMax

    % Plot the final B-spline curve
    plot(curvePoints(1, :), curvePoints(2, :), 'DisplayName', sprintf('B-spline Curve (k=%d)', k), ...
        'LineWidth', 1.5)
    hold on
end

% Function to generate the knot vector
function kv = generateKnotVector(n, k)
    % Generates a uniform knot vector
    % Inputs:
    % n - number of control points
    % k - order of the B-spline
    % Output:
    % kv - the knot vector
    
    kv = zeros(1, n+k+1);
    kv((k+1):(n+1)) = (k+1 : n+1) - k;
    kv((n+2):(n+k+1)) = n - k + 2;
end

% Function to calculate the basis function N(i,k) at a given parameter t
function basisValue = basisFunction(i, k, t, kv)
    % Inputs:
    % i - index of the basis function
    % k - order of the basis function
    % t - parameter value
    % kv - knot vector
    % Output:
    % basisValue - value of the basis function N(i,k) at t
    
    if (k == 1)
        % Base case: piecewise constant basis function
        if (kv(i) <= t && t < kv(i+1))
            basisValue = 1;
        else
            basisValue = 0;
        end
    else
        % Recursive definition of B-spline basis functions
        d1 = kv(i+k-1) - kv(i);
        d2 = kv(i+k) - kv(i+1);
        
        s1 = 0; s2 = 0;
        
        if (d1 ~= 0)
            s1 = ((t - kv(i)) * basisFunction(i, k-1, t, kv)) / d1;
        end
        
        if (d2 ~= 0)
            s2 = ((kv(i+k) - t) * basisFunction(i+1, k-1, t, kv)) / d2;
        end
        
        basisValue = s1 + s2;
    end
end

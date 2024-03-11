function [x, t, u] = calctemp(tmax, nt, xmax, nx, method, timeData, tempData)
% Function for modelling temperature in a space shuttle tile
% D N Johnston  14/02/24
%
% Input arguments:
% tmax   - maximum time (s)
% nt     - number of timesteps
% xmax   - total thickness (m)
% nx     - number of spatial steps
% method - solution method ('forward', 'backward' etc)
% timeData - time vector for surface temperatures (s)
% tempData - surface temperature vector (C or K)
%
% Return arguments:
% x      - distance vector (m)
% t      - time vector (s)
% u      - temperature matrix (C or K)
%
% For example, to perform a  simulation with 501 time steps
%   [x, t, u] = calctemp(4000, 501, 0.05, 21, 'forward', timeData, tempData);
%

% Set material properties and derived values (LI-900)
% Obtained from NASA document: Structures and Materials: Space Shuttle Tiles, Grades 5-8 - NASA
% Note that we're assuming constant properties.
thermCon = 0.0484; % W/m K; 0.028 BTU/ft/hr/F at 70F and 1 atm
density  = 144;    % kg/m^3; 9 lb/ft^3
specHeat = 628;    % J/kg/K; 0.15 Btu/lb/F

% Initialise everything.
dt = tmax / (nt-1);
t = (0:nt-1) * dt;
dx = xmax / (nx-1);
x = (0:nx-1) * dx;
u = zeros(nt, nx);

% Use interpolation to get outside temperature at time vector t 
% and store it as left-hand boundary vector L.
L = interp1(timeData, tempData, t, "linear", "extrap");

% set initial conditions equal to boundary temperature at t=0.
u(1, :) = L(1);

% Select method and run simulation.
switch method
    case 'forward'
        u(:, 1) = L; % Outside boundary condition
        u(:, nx) = 0; % Inside boundary condition; set to zero as a starting point.
        % You need to put your solution code here.
        p = thermCon * dt / dx^2;
        for n = 1:nt-1
            % calculate internal values using forward differencing
            for i = 2:nx-1
                u(n+1,i) = (1 - 2 * p) * u(n,i) + p * (u(n,i-1) + u(n,i+1));
            end
        end

    case 'dufort-frankel'
        u(:, 1) = L;
        u(:, nx) = u(:, 1);
        p = thermCon * dt / dx^2;
        for n = 2:nt-1
            for i = 2:nx-1
                u(n+1,i) = ((2 * p) * (u(n,i-1) + u(n,i+1)) + (1 - 2 * p) * u(n-1,i)) / (1 + 2 * p);
            end
        end

    case 'backward'
        u(:, 1) = L;
        A = diag((1 + 2*p)*ones(1, nx-2)) + diag(-p*ones(1, nx-3), 1) + diag(-p*ones(1, nx-3), -1);
        A(1,2) = -2*p;
        A(end,end-1) = -2*p;
        for n = 1:nt-1
            B = u(n, 2:nx-1)';
            B(1) = B(1) + p*u(n+1, 1);
            B(end) = B(end) + p*u(n+1, nx);
            u(n+1, 2:nx-1) = A\B; 
        end
        
    case 'crank-nicolson'


    otherwise
        error (['Undefined method: ' method])
end



    
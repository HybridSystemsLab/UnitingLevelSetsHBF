function xdot = fU(x)
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: f.m
%--------------------------------------------------------------------------
% Project: Simulation of the Heavy Ball method for finding the nearest
% non-unique minimum. This is non-hybrid, for now.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00
   
% The global variables
global lambda_0 gamma_0 lambda_1 gamma_1 

% state
z1 = x(1:100,1);
z2 = x(101:200,1);
q = x(201,1);

% Black box: L could be anything. 
y1 = GradientL(z1);

if (q == 0) 
    u = - lambda_0*z2 - gamma_0*y1; 
elseif (q == 1)
    u = - lambda_1*z2 - gamma_1*y1; 
end

xdot = [z2;u;0]; 
end
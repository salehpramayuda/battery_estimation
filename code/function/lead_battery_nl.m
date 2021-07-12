function [dx, y] = lead_battery_nl(t, x, u, a1, a2, a3, Cm1, Rm1, ...
R0dm1, R0cm1, varargin)
%% Setup
Q_e = 2.9*3600;
if u>=0
    R0m1 = R0cm1;
else
    R0m1 = R0dm1;
end

%% Output
y = a1*x(1)^2+a2*x(1)+a3 + x(2) + u/R0m1;

%% Derivative
dx(1) = u/Q_e;
dx(2) = Cm1*(u-x(2)*Rm1);


end
% validation for the basic EKF implementation
clc; clear
addpath('./../modell_validation/')
%% Simulation Parameter
% biphasic current
t_sim = 1000;
Delta = 0.04;
[t, i] = make_bipulse(0.5, t_sim, 20, 15, 10, 30, 100, Delta, 0.25);
cur_input = timeseries(i,t);
if(0)
    figure(); grid on
    plot(t, i); xlabel('Time [s]'); ylabel('Current [A]');
    legend('i_{sim}');
end

% Battery parameter
Rp      = 1.8e3;
R0      = 0.75e3;
Cp      = 300e-3;
Cg      = 2*3.3e-4;
Qe      = 2.9*3600;
soc0    = 0.3;

% calculate SoC parameter
u_oc    = 2*[9.75, 10.1, 10.3, 10.55, 10.75, 10.9, 11.2, 11.3, 11.45, ...
    11.55, 11.7, 11.8, 12.05, 12.15, 12.25, 12.35, 12.46, 12.55, ...
    12.60, 12.60, 12.65]';
soc     = 0:0.05:1;
phi     = [soc'.^2 soc' ones(21,1)];
a       = (phi'*phi)\phi'*u_oc;
clear u_oc soc phi

%% make voltage data
filename = 'sim_battery_wo_cap';
systemModel = load_system(filename);
modelConfig = getActiveConfigSet(systemModel).copy;
simResult = sim(filename, modelConfig);
t = simResult.tout; simResult = simResult.simout;

% Werte aus der Simulation
u_L = simResult(:,1)';
i_b = simResult(:,2)';
u_rc = simResult(:,3)';
soc = simResult(:,4)';

time_est    = 400; [~, ind_est] = min(abs(t-time_est));
cur_est     = timeseries(i_b(1:ind_est), t(1:ind_est));
volt_est    = timeseries(u_L(1:ind_est), t(1:ind_est));

%% setup EKF
% Initial state vector
x0 = [1/Cp; 1/Rp; 1/R0; soc0; 0];
x0(5) = volt_est.Data(1)-cur_est.Data(1)/x0(3)-a(1)*x0(4)^2-a(2)*x0(4)-a(3);

% Rough estimation for initial covariance-Matr.
P0 = eye(5); P0(4,4) = soc0^2; P0(5,5) = (x0(5)*0.3)^2;
for j = 1:3
    P0(j,j) = (x0(j)*0.3)^2;
end
% Variance of system noise and measurement noise
W_k = zeros(length(x0));
Z_k = 0;
clear j

%% run EKF
systemModel = load_system('parameter_estimation');
modelConfig = getActiveConfigSet(systemModel).copy;
simResult = sim('parameter_estimation', modelConfig);

t_sim = simResult.tout;
simResult = reshape(simResult.simout,7,[]);
x_conv = simResult(1:5,end); 







% validation for the basic EKF implementation
clc; clear
addpath('./../modell_validation/')
%% Simulation Parameter
% biphasic current
t_sim = 1000;
Delta = 0.04;
[t, i] = make_bipulse(0.5, t_sim, 20, 15, 10, 30, 100, Delta, 0.25);
cur_input = timeseries(i,t);

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
filename = 'sim_battery';
systemModel = load_system(filename);
modelConfig = getActiveConfigSet(systemModel).copy;
simResult = sim(filename, modelConfig);
t = simResult.tout; simResult = simResult.simout;

% Werte aus der Simulation
u_L = simResult(:,1)';
i_b = simResult(:,2)';
u_rc = simResult(:,3)';
soc = simResult(:,4)';
i_L = simResult(:,5)';
i_C = simResult(:,6)';

time_est    = 400; [~, ind_est] = min(abs(t-time_est));
cur_est     = timeseries(i_b(1:ind_est), t(1:ind_est));
volt_est    = timeseries(u_L(1:ind_est), t(1:ind_est));
volt_est_kp1= timeseries(u_L(2:ind_est+1), t(1:ind_est));

%% setup EKF
Delta = cur_est.Time(2)-cur_est.Time(1);
Z_k = 0;

% Initial state vector
x0 = [2/Cp; 1/Rp; 2/R0; 0.30; 0; volt_est.Data(1)];

P0 = eye(6); P0(4,4) = x0(4)^2; P0(5,5) = (x0(5)*0.3)^2; P0(6,6) = Z_k;
for i = 1:3
    P0(i,i) = (x0(i)*0.3)^2;
end

W_k = zeros(6);

clear i

%% Simulate EKF

systemModel = load_system('parameter_estimate');
modelConfig = getActiveConfigSet(systemModel).copy;
simResult = sim('parameter_estimate', modelConfig);

t_sim = simResult.tout;
simResult = simResult.simout;
x_conv = simResult(end, :);






clc;clear;
%% Setup
addpath('./../..')
load('../../src/210609_data');
addpath('./function/')
main_data = data{11};
volt_ind = find(ismember(getElementNames(main_data),'Voltage(V)'));
if isempty(volt_ind)
    disp('Dataset does not contain measured voltage');
    return
end
cur_ind = find(ismember(getElementNames(main_data),'Current(A)'));
if isempty(cur_ind)
    disp('Dataset does not contain measured current');
    return
end

voltage = main_data{volt_ind}.Values;
current = main_data{cur_ind}.Values;
time_mes = voltage.Time(end);

% Split measured data into estimation and validation
time_est = 350;
time_val = time_mes - time_est;

[~,ind_val] = min(abs(voltage.Time-time_est));

volt_est = timeseries(voltage.Data(1:ind_val),voltage.Time(1:ind_val));
volt_val = timeseries(voltage.Data(ind_val+1:end),...
    voltage.Time(ind_val+1:end)-voltage.Time(ind_val+1));

cur_est = timeseries(current.Data(1:ind_val),current.Time(1:ind_val));
cur_val = timeseries(current.Data(ind_val+1:end),...
    current.Time(ind_val+1:end)-current.Time(ind_val+1));


%% Find Variance of zeta Z_k
% for i as input
[m,ind] = min(abs(voltage.Time-55));
static_volt = voltage.Data(1:ind);
Z_b_k = var(static_volt);

% % for u as input
[m,ind] = min(abs(current.Time-55));
static_curr = current.Data(1:ind);
Z_k = var(static_curr);
clear m

%% SOC-Graph
u_oc = 2*[9.75, 10.1, 10.3, 10.55, 10.75, 10.9, 11.2, 11.3, 11.45, ...
    11.55, 11.7, 11.8, 12.05, 12.15, 12.25, 12.35, 12.46, 12.55, ...
    12.60, 12.60, 12.65]';
soc = 0:0.05:1;
phi = [soc'.^2 soc' ones(21,1)];
a = inv(phi'*phi)*phi'*u_oc;
% figure(1);
% plot(soc, u_oc); hold on
% plot(soc, phi*param);
% grid on; hold off

%exportgraphics(gcf, 'graph/soc_kurve.pdf', 'ContentType', 'Vector');
%% first EKF
Delta = 8e-3; 
Qe = 2.9*3600;
u_oc0 = a(1)*0.5^2+a(2)*0.5+a(3);
Cg = 2*3.3e-4; % taken from Vrontos-Paper

% setup state-ekf
% x = [1/C; 1/R; 1/R0; SOC; U_rc; i_l, i_b]
% Initial value: x1_0 = 2/790.27; x2_0 = 1/0.0252; x3_0 = 1/0.0172;
% x4_0 = 50%; x5_0 = 0; x6_0 = 0; x7_0 = i_b(0)

x_0 = [2/790.27; 1/0.0252; 1/0.0172; 0.5; 0; 0; cur_est.Data(1)];

% set covariance matrix
P_0 = eye(6);
for i = 1:3
    P_0(i,i) = (0.3*x_0(i))^2;
end
P_0(4,4) = 0.5^2; P_0(5,5) = (volt_est.Data(1)*0.05)^2;
P_0(6,6) = (x_0(7)+ (volt_est.Data(2)-volt_est.Data(1))*Cg)^2;
P_0(7,7) = (Z_k)^2;

% error variance matrix
W_k = zeros(7,7); W_k(6,6) = 0.1;

%% Setup for last EKF

% last-ekf
% x_b = Q_cg with initial value Q_cg = u_L(0)*Cg
x_b_0 = volt_est.Data(1)*Cg;

% covariance matrix for last-ekf
P_b_0 = (0.3*x_b_0)^2;

% error variance matrix for last-ekf
W_b_k = 0;
%% Simulate EKF
filename = 'akku_schaetzung';
systemModel = load_system(filename);
modelConfigSet = getActiveConfigSet(systemModel);
modelConfig = modelConfigSet.copy;
simResult = sim(filename, modelConfig);

x = simResult.simout(1:8,:);
x_conv = x(:,end);

%% Validate Model

filename = 'validierung';
systemModel = load_system(filename);
modelConfigSet = getActiveConfigSet(systemModel);
modelConfig = modelConfigSet.copy;
simResult = sim(filename, modelConfig);
%%
t_val = simResult.tout(5:end);
y_hat = simResult.simout(5:end,1);
error = simResult.simout(5:end,3);

figure();
plot(t_val, y_hat); hold on; grid on
%plot(volt_val);
plot(cur_val); % for u = u_l, y = i_bat
hold off
legend('simulated', 'measured');

figure();
plot(t_val, error);grid on;

fit = 100*(1-goodnessOfFit(y_hat, simResult.simout(5:end,4), 'NRMSE'))

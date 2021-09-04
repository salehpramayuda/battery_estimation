clc;clear;
%% Setup
addpath('./../..')
load('../../src/210609_data');
addpath('./../../function/')
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
time_est_1 = 140;
time_est_2 = 400; time_dekf = time_est_2-time_est_1;
time_val = time_mes - time_est_2;

[~,ind_1] = min(abs(voltage.Time-time_est_1));
[~,ind_2] = min(abs(voltage.Time-time_est_2));

volt_est_first = timeseries(voltage.Data(1:ind_1),voltage.Time(1:ind_1));
volt_est_sec = timeseries(voltage.Data(ind_1+1:ind_2),voltage.Time(ind_1+1:ind_2)-...
    voltage.Time(ind_1+1));
volt_val = timeseries(voltage.Data(ind_2+1:end),...
    voltage.Time(ind_2+1:end)-voltage.Time(ind_2+1));

cur_est_first = timeseries(current.Data(1:ind_1),current.Time(1:ind_1));
cur_est_sec = timeseries(current.Data(ind_1+1:ind_2),current.Time(ind_1+1:ind_2)-...
    current.Time(ind_1+1));
cur_val = timeseries(current.Data(ind_2+1:end),...
    current.Time(ind_2+1:end)-current.Time(ind_2+1));

clear voltage current volt_ind cur_ind ind_1 ind_2 main_data
%% Find Variance of zeta Z_k
% for i as input
[~,ind] = min(abs(volt_est_first.Time-55));
static_volt = volt_est_first.Data(1:ind);
Z_b_k = std(static_volt)^2;

% % for u as input
[~,ind] = min(abs(cur_est_first.Time-55));
static_curr = cur_est_first.Data(1:ind);
Z_k = std(static_curr)^2;

clear ind static_curr static_volt
%% SOC-Graph
u_oc = 2*[9.75, 10.1, 10.3, 10.55, 10.75, 10.9, 11.2, 11.3, 11.45, ...
    11.55, 11.7, 11.8, 12.05, 12.15, 12.25, 12.35, 12.46, 12.55, ...
    12.60, 12.60, 12.65]';
soc = 0:0.05:1;
phi = [soc'.^2 soc' ones(21,1)];
a = (phi'*phi)\phi'*u_oc;

clear soc phi
%% EKF first phase
Delta = 8e-3; 
Qe = 2.9*3600;
u_oc0 = a(1)*0.5^2+a(2)*0.5+a(3);

% setup state-ekf
% x = [1/C; 1/R; 1/R0; SOC; U_rc]

x0 = [2/790.27; 1/0.0252; 1/0.0172; 0.5; 0];
P0 = eye(5); P0(1,1)=(0.03*x0(1))^2;P0(2,2)=(0.03*x0(2))^2;P0(3,3)=(0.03*x0(3))^2;
P0(4,4)=0.5^2;P0(5,5)=(volt_est_first.Data(1)-u_oc0-cur_est_first.Data(1)/x0(3))^2;
W_k = zeros(5,5);

%% Simulate EKF
systemModel = load_system('dekf_estimate_attempt_3_early');
modelConfig = getActiveConfigSet(systemModel).copy;
simResult = sim('dekf_estimate_attempt_3_early', modelConfig);

x = simResult.simout(:,1:5);
x_conv = x(end,:);

%% (D-)EKF second phase
Cg = 2*3.3e-4; % taken from Vrontos-Paper

% setup state-ekf with x = [1/C; 1/R; 1/R0; SOC; U_rc; i_l, i_b]
i_l0 = (volt_est_sec.Data(2)-volt_est_sec.Data(1))/Delta*Cg+cur_est_sec.Data(1);
x0 = [x_conv'; i_l0; cur_est_sec.Data(1)];

P0 = eye(7);
for i = 1:3
    P0(i,i) = (0.3*x0(i))^2;
end
P0(4,4) = 0.5^2; P0(5,5) = (volt_est_sec.Data(1)*0.05)^2;
P0(6,6) = (x0(6)*0.1)^2;
P0(7,7) = (Z_k)^2;
W_k = zeros(7,7); W_k(6,6) = 0.1;

% second-ekf
% x_b = Q_cg with initial value Q_cg = u_L(0)*Cg

x_b_0 = volt_est_sec.Data(1)*Cg;
P_b_0 = (0.3*x_b_0)^2;
W_b_k = 0;

%% Simulate EKF

systemModel = load_system('dekf_estimate_attempt_3');
modelConfig = getActiveConfigSet(systemModel).copy;
simResult = sim('dekf_estimate_attempt_3', modelConfig);

x = simResult.simout;
x_conv_2 = x(:,end)';


%% Validate Model
i_l_ = (volt_val.Data(2)-volt_val.Data(1))/Delta*Cg+cur_val.Data(1);

x0 = [x_conv_2(6); cur_val.Data(1)];
P0 = eye(2); P0(2,2) = Z_k^2; P0(1,1) = (x0(1)-i_l_)^2;
W_k = zeros(2,2); W_k(1,1) = 0.1;
%%
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

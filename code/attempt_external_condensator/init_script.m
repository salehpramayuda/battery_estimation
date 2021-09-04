clc;clear;
%% Load and Extract Data
load('./../src/210609_data');
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
time_est = 400;
time_val = time_mes - time_est;

[~,ind_val] = min(abs(voltage.Time-time_est));

volt_est = timeseries(voltage.Data(1:ind_val),voltage.Time(1:ind_val));
volt_val = timeseries(voltage.Data(ind_val+1:end),...
    voltage.Time(ind_val+1:end)-voltage.Time(ind_val+1));

cur_est = timeseries(current.Data(1:ind_val),current.Time(1:ind_val));
cur_val = timeseries(current.Data(ind_val+1:end),...
    current.Time(ind_val+1:end)-current.Time(ind_val+1));


clear data main_data volt_ind cur_ind voltage current time_mes ind_val 
%% Find Variance of Zeta Z_k
% for i as input
t_static = 55;  % time before stimulations
[~,ind] = min(abs(volt_est.Time - t_static));
static_volt = volt_est.Data(1:ind);
Z_k = var(static_volt);

clear t_static ind static_volt
%% SOC-Graph
u_oc = 2*[9.75, 10.1, 10.3, 10.55, 10.75, 10.9, 11.2, 11.3, 11.45, ...
    11.55, 11.7, 11.8, 12.05, 12.15, 12.25, 12.35, 12.46, 12.55, ...
    12.60, 12.60, 12.65]';
soc = 0:0.05:1;
phi = [soc'.^2 soc' ones(21,1)];
a = (phi'*phi)\phi'*u_oc;

figure()
plot(soc, u_oc, 'r'); hold on
plot(soc, phi*a, 'b'); hold off
grid on
legend('Measured', 'Estimated');
ylabel('Open Circuit Voltage U_{OC} [V]');
xlabel('State of Charge');

clear u_oc soc phi
return
%% Initial Values
% Important parameter
Delta = cur_est.Time(2)-cur_est.Time(1);
Qe = 2.9*3600;

% Initial state vector
x0 = [2/(2.9*3600/12); 1/42e-3/2; 1/42e-3/2; 0.65; 0];
x0(5) = volt_est.Data(1)-cur_est.Data(1)/x0(3)-a(1)*x0(4)^2-a(2)*x0(4)-a(3);

% failed initial state values:
% u_rc,0 = 0 -> should be somewhat calculated and not 0
% soc,0 works well with 0.5 or 0.75

% Rough estimation for initial covariance-Matr.
P0 = eye(5); P0(4,4) = 0.65^2; P0(5,5) = (x0(5)*0.3)^2;
for i = 1:3
    P0(i,i) = (x0(i)*0.3)^2;
end

% Variance of system noise
W_k = zeros(length(x0));

clear i

%% Simulate EKF to Extract Converged Parameter
systemModel = load_system('parameter_estimation');
modelConfig = getActiveConfigSet(systemModel).copy;
simResult = sim('parameter_estimation', modelConfig);

t_sim = simResult.tout;
simResult = reshape(simResult.simout,7,[]);
x_conv = simResult(1:5,end); 
%%
if(1)
    param = simResult(1:5,:,:);
    y_hat = simResult(6,:,:);
    error = simResult(7,:,:);
    
    figure(1);
    subplot(3,2,1); plot(t_sim, 1/param(1,:)'); ylabel('Capacitor C'); grid on
    subplot(3,2,2); plot(t_sim, 1/param(2,:)'); ylabel('Resistance R');grid on
    subplot(3,2,3); plot(t_sim, 1/param(3,:)'); ylabel('Resistance R0');grid on
    subplot(3,2,4); plot(t_sim, param(5,:)'); ylabel('RC-Voltage'); grid on
    subplot(3,2,5:6); plot(t_sim, param(4,:)'); ylabel('State of Charge');grid on
    title('Estimation Results');
    
    figure(2);
    subplot(2,1,1);
    plot(t_sim, y_hat); hold on; plot(volt_est.Time, volt_est.Data);
    ylabel('Voltage [V]'); xlabel('Time [s]');
    legend('Estimated', 'Measured'); hold off; grid on
    subplot(2,1,2); 
    plot(t_sim, error);
    ylabel('Error Voltage [V]'); xlabel('Time [s]'); grid on
    
    clear y_hat error
end

%% Simulate Battery Model for Validation
Cg = 2*3.3e-4; % taken from Vrontos-Paper
x_conv = [volt_val.Data(1); x_conv];

% systemModel = load_system('simulate_battery');
% modelConfig = getActiveConfigSet(systemModel).copy;
% simResult = sim('simulate_battery', modelConfig);



% Simulation script
clear; clc;
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

if(0)
    figure()
    plot(soc, u_oc, 'r'); hold on
    plot(soc, phi*a, 'b'); hold off
    grid on
    legend('Measured', 'Estimated');
    ylabel('Open Circuit Voltage U_{OC} [V]');
    xlabel('State of Charge');
end
clear u_oc soc phi

%% Run simulation
w_cap = true;
if(w_cap)
    filename = 'sim_battery';
else
    filename = 'sim_battery_wo_cap';
end
systemModel = load_system(filename);
modelConfig = getActiveConfigSet(systemModel).copy;
simResult = sim(filename, modelConfig);
t = simResult.tout; simResult = simResult.simout;

% Werte aus der Simulation
u_L = simResult(:,1)';
i_b = simResult(:,2)';
u_rc = simResult(:,3)';
soc = simResult(:,4)';
if(w_cap)
    i_L = simResult(:,5)';
    i_C = simResult(:,6)';
end

%% Run with capacitor
if(w_cap)
    x_conv = [1/Cp, 1/Rp, 1/R0, soc0, 0, 0, a(1)*soc0^2+a(2)*soc0+a(3)]';
    cur_val = cur_input;
    systemModel = load_system('validate_model_w_cap');
    modelConfig = getActiveConfigSet(systemModel).copy;
    simValidate = sim('validate_model_w_cap', modelConfig);
else
    x_conv = [1/Cp, 1/Rp, 1/R0, soc0, 0]';
    cur_val = cur_input;
    systemModel = load_system('validate_model_wo_cap');
    modelConfig = getActiveConfigSet(systemModel).copy;
    simValidate = sim('validate_model_wo_cap', modelConfig);
end
t_val = simValidate.tout;
u_L_sim = simValidate.simout(:,1)';
soc_sim = simValidate.simout(:,2)';
u_rc_sim = simValidate.simout(:,3)';

figure();
subplot(2,1,1); hold on
plot(t, u_L, 'g');
plot(t_val, u_L_sim, 'r--'); hold off
grid on; legend('circuit model', 'simulated model');
subplot(2,1,2); plot(t, i_b, 'b'); grid on;

figure();
subplot(2,1,1); hold on
plot(t, soc, 'g');
plot(t_val, soc_sim, 'r--'); hold off
grid on; legend('circuit model', 'simulated model');
subplot(2,1,2); hold on
plot(t, u_rc, 'b');
plot(t_val, u_rc_sim, 'r--'); hold off
grid on; legend('circuit model', 'simulated model');

if w_cap
    i_b_sim = simValidate.simout(:,4)';
    i_c_sim = simValidate.simout(:, 5)';
    figure(); subplot(2,1,1);hold on
    plot(t, i_C, 'g'); plot(t_val, i_c_sim, 'r--'); hold off
    grid on; legend('circuit model', 'simulated model');
   
    subplot(2,1,2); hold on
    plot(t, i_b, 'b'); plot(t_val, i_b_sim, 'r--'); hold off
    grid on; legend('circuit model', 'simulated model');
end


%% Load and Run Simulink
clc
filename = 'akku_schaetzung';
systemModel = load_system(filename);
modelConfigSet = getActiveConfigSet(systemModel);
modelConfig = modelConfigSet.copy;
simResult = sim(filename, modelConfig);
%%
%x = simResult.simout(:,1:5);
x = simResult.simout_check(:,3:7);
u = simResult.simout(:,6);
y = simResult.simout(:,7);
t = simResult.tout;
y_hat = simResult.simout_check(:,1);
y_hat_ = simResult.simout_check(:,2);

%% Plot
g = figure('Name', 'EKF Zustände');
g.NextPlot = 'add';
ax = axes;
ht = title('State = (1/C, 1/R, 1/R0, SoC, U_RC)');
ax.Visible ='off';
ht.Visible = 'on';
for i = 1:5
    subplot(3,2,i);
    plot(t, x(:,i),'r', 'LineWidth', 1.5);
    xlabel('Time [s]');
    grid();
end
subplot(3,2,1);ylabel('x1: 1/C [1/F]');
subplot(3,2,2);ylabel('x2: 1/R [1/\Omega]');
subplot(3,2,3);ylabel('x3: 1/R0 [1/\Omega]');
subplot(3,2,4);ylabel('x4: SoC');
subplot(3,2,5);ylim([-0.001,0.02]);ylabel('x5: U_{RC} [V]');
subplot(3,2,6);
plot(t, u, 'g', 'LineWidth', 1.5);
grid();ylabel('Current [A]');xlabel('Time [s]');

%param_conv = [1/x(end,1), 1/x(end,2), 1/x(end,3)];
%save('converged_param.mat', 'param_conv');
%exportgraphics(gcf, 'ekf_result.pdf', 'ContentType', 'Vector');

% %% 2nd Plot
% g = figure('Name', 'Voltages');
% g.NextPlot = 'add'; ax = axes;
% ht = title('Voltages on Battery'); ax.Visible='off';
% ht.Visible='on';
% subplot(2,2,1);plot(t, y, 'r');title('Terminal Voltage');
% ylabel('U_L [V]');xlabel('Time [s]');grid on;
% 
% subplot(2,2,2);plot(t, [x(:,4).^2 x(:,4) ones(length(x(:,4)),1)]*param, 'g');
% title('Open Circuit Voltage');
% ylabel('U_{OC} [V]');xlabel('Time [s]');grid on;
% 
% subplot(2,2,3);plot(t, x(:,5), 'b');title('Voltage over C');
% ylabel('U_{RC} [V]');xlabel('Time [s]');grid on;
% 
% subplot(2,2,4);plot(t, u./x(:,3), 'r');title('Voltage over R_0');
% ylabel('U_{R0} [V]');xlabel('Time [s]');grid on;
% 
% exportgraphics(gcf, 'battery_voltages.pdf', 'ContentType', 'Vector');


%% Check output
fit = 100*(1-goodnessOfFit(y_hat, y, 'NRMSE'));
h = figure('Name', 'Check Output');
grid on; hold on;
plot(t, y, 'g', 'LineWidth', 1.5);
plot(t, y_hat, 'r--');
legend({'y_{mess}','y_{sim param-const}'}); title('Check Output');
xlabel('Time [s]');ylabel('Voltage [V]');

disp(fit);

%%
f= figure('Name', 'Output Error');
plot(t, y-y_hat);
grid on
xlabel('Time [s]'); ylabel('Voltage [V]');











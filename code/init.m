clc; clear;
%% Setup
load('Parameters_and_Signals.mat');
main_data = data{5};
volt_ind = find(ismember(getElementNames(main_data),'Voltage (V)'));
if isempty(volt_ind)
    disp('Dataset does not contain measured voltage');
    return
end
cur_ind = find(ismember(getElementNames(main_data),'Current (A)'));
if isempty(cur_ind)
    disp('Dataset does not contain measured current');
    return
end

voltage = main_data{volt_ind}.Values;
current = main_data{cur_ind}.Values;

%% Clean Data
addpath('./function')
[volt_clean, place] = clean(voltage,8);
cur_clean = timeseries(current.Data(place), current.Time(place));

%% Plot
figure('Name', 'Voltage-Current-Dataset 5');grid on;
yyaxis left
plot(volt_clean,'Color',[0.5,0.15,0.45]);
ylabel('Voltage [V]');
xlabel('Time [s]');
% yyaxis right
% plot(cur_clean.Time,cur_clean.Data, 'b');
% ylabel('Current [A]');
title('Dataset 5');

%% Init for Simulink
Delta = volt_clean.Time(2)-volt_clean.Time(1);
% x1 = Spannung über Kondensator; x2 = Klemmenspannung
% x3 = 1/C; x4 = 1/R (Elektrochemischer Widerstand)
% Initial state x[0|-1]
% Uc = 0 (not yet charged), UL = y,
% 1/C = U/(I*t), 1/R = 1/42e-3 (Laut Datenblatt)
% x0 = [0, volt_clean.Data(1), 12/(2.9*3600), 1/42e-3]';
% P0 = eye(4); P0(1,1) = (x0(2)*exp(-Delta*x0(3)*x0(4)))^2; P0(2,2) = 0.02^2;
% P0(3,3) = (0.03*x0(3))^2; P0(4,4) = (0.03*x0(4))^2;
% W_k = 0*eye(5);
% Z_k = 0.02;

% Corrected EKF Implementation
% x1 = 1/C; x2 = 1/R; x3 = 1/R0; x4 = Uoc
% x5 = U (Spannung über Kondensator)
% Initial state x[0|-1]
% 1/C = U/(I*t), 1/R = 1/42e-3 (Datenblatt), 
% 1/R0 = y_k/(x4_k-x5_k-u)
% Uoc = 12V (Datenblatt),Uc = 0 (not yet charged)

% x0 = [12/(2.9*3600), 1/42e-3, cur_clean.Data(1)/(25.5-0-volt_clean.Data(1)),...
%     25, 0]';
% x0(1) = 1/790.27; x0(2) = 1/0.0126;
% 
% P0 = eye(5); P0(1,1)=(0.03*x0(1))^2;P0(2,2)=(0.03*x0(2))^2;P0(3,3)=(0.03*x0(3))^2;
% P0(4,4)= (0.3*x0(4))^2;
% P0(5,5)=(x0(4)-volt_clean.Data(1)-cur_clean.Data(1)/x0(3))^2;

%% 3. Version of EKF
% x1 = 1/C; x2 = 1/R; x3 = 1/R0; x4 = 1/Qe;
% x5 = SOC; x6 = U (Voltage over Capacitor)
% Initial state
% x1_0 = 1/790.27; x2_0 = 1/0.0126; x3_0 = 1/0.0086;
% x4_0 = 1/(2.9*3600); x5_0 = 50%; x6_0 = 0;

x0 = [1/790.27; 1/0.0126; 1/0.0086; 1/(2.9*3600); 0.5; 0];
alpha = 5/2; beta = 41/4;

P0 = eye(5); P0(1,1)=(0.03*x0(1))^2;P0(2,2)=(0.03*x0(2))^2;P0(3,3)=(0.03*x0(3))^2;
P0(4,4)=(x0(4)*(1-100/97))^2;P0(5,5)=0.5^2;P0(6,6)=(volt_clean.Data(1)*0.05)^2;

W_k = 0.0*eye(6);
Z_k = 0.02;

%% Plot
% filename = 'akku_schaetzung';
% systemModel = load_system(filename);
% modelConfigSet = getActiveConfigSet(systemModel);
% modelConfig = modelConfigSet.copy;
% simResult = sim(filename, modelConfig);
% 
% 
% h = figure('Name', 'Extended Kalman Filter');
% h.NextPlot = 'add';
% ax = axes;
% ht = title('Zustand = (1/C, 1/R, 1/R0, Uoc, U)');
% ax.Visible ='off';
% ht.Visible = 'on';
% for i = 1:5
%     subplot(3,2,i);
%     plot(simResult.tout, simResult.simout(:,i));
%     grid();
% end
% subplot(3,2,1);ylabel('x1: 1/C');
% subplot(3,2,2);ylabel('x1: 1/R');
% subplot(3,2,3);ylabel('x1: 1/R0');
% subplot(3,2,4);ylabel('x1: Uoc');
% subplot(3,2,5);ylabel('x1: U');
% subplot(3,2,6);
% plot(cur_clean);
% grid();
% ylabel('Strom');
% 
% exportgraphics(gcf, 'zustand_yasim.pdf', 'ContentType', 'Vector');

% function x_k  = extended_kf(u_k,y_k,x0, P0, W_k, Z_k, Delta)
%     % Extended Kalman Filter for estimating Capacity and U_L
%     
%     % Initialisierung Prädiktion bei k = 0
%     persistent P_k_p x_k_p
%     if isempty(P_k_p)
%         P_k_p = P0;
%     end
%     if isempty(x_k_p)
%         x_k_p = x0;
%     end
%     
%     % y_hat_k, C, K, x_k, P_k
%     y_hat_k = x_k_p(3)*(x_k_p(4)-x_k_p(5)-u_k);
%     C_k = [0, 0, x_k_p(4)-x_k_p(5)-u_k, x_k_p(3), -x_k_p(3)];
%     K_k = P_k_p*C_k'/(Z_k + C_k*P_k_p*C_k');
%     x_k = x_k_p + K_k*(y_k - y_hat_k);
%     P_k = P_k_p - K_k*C_k*P_k_p;
%     
%     % Calculate A, prediction of x and P
%     A_k = [[eye(4), zeros(4,1)]; [Delta*(x_k(3)*(x_k(4)-u_k)-x_k(5)*...
%         (x_k(2)+x_k(3))), -Delta*x_k(1)*x_k(5), Delta*x_k(1)*...
%         (x_k(4)-u_k-x_k(5)), Delta*x_k(1)*x_k(3),1-Delta*x_k(1)*(x_k(2)+x_k(3))]]; 
%     x_kp1_p = [x_k(1);x_k(2);x_k(3);x_k(4);x_k(5)+Delta*x_k(1)*...
%         (x_k(3)*(x_k(4)-u_k)-x_k(5)*(x_k(2)+x_k(3)))];
% 
%     P_kp1_p = A_k*P_k*A_k' + W_k;
%     
%     % Update persistent variable
%     x_k_p = x_kp1_p;
%     P_k_p = P_kp1_p;
% end

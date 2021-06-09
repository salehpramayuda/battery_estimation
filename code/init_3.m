%% Setup
load('210520_Versuch');
main_data = data{1};
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

% %% Clean Data (not needed, data is already clean)
% [volt_clean, place] = clean(voltage,8);
% cur_clean = timeseries(current.Data(place), current.Time(place));


%% Find Variance of zeta Z_k
[m,ind] = min(abs(voltage.Time-100));
static_volt = voltage.Data(1:ind);
Z_k = var(static_volt);

%% SOC-Kennlinie 
u_oc = 2*[9.75, 10.1, 10.3, 10.55, 10.75, 10.9, 11.2, 11.3, 11.45, ...
    11.55, 11.7, 11.8, 12.05, 12.15, 12.25, 12.35, 12.46, 12.55, ...
    12.60, 12.60, 12.65]';
soc = 0:0.05:1;
phi = [soc'.^2 soc' ones(21,1)];
param = inv(phi'*phi)*phi'*u_oc;
%% EKF (with Q_0 as parameter/state)
% x1 = 1/C; x2 = 1/R; x3 = 1/R0; x4 = 1/Q_0;
% x5 = SOC; x6 = U (Voltage over Capacitor)
% Initial state
% x1_0 = 1/790.27; x2_0 = 1/0.0126; x3_0 = 1/0.0086;
% x4_0 = 1/(2.9*3600); x5_0 = 50%; x6_0 = 0;

x0 = [2/790.27; 1/2/0.0126; 1/2/0.0086; 1/(2.9*3600); 0.5; 0];
a1= param(1); a2=param(2); a3=param(3);

P0 = eye(6); P0(1,1)=(0.03*x0(1))^2;P0(2,2)=(0.03*x0(2))^2;P0(3,3)=(0.03*x0(3))^2;
P0(4,4)=(x0(4)*(1-100/97))^2;P0(5,5)=0.5^2;P0(6,6)=(volt_clean.Data(1)*0.05)^2;

W_k = 0.0*eye(6);

%%
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

%exportgraphics(gcf, 'zustand_yasim.pdf', 'ContentType', 'Vector');
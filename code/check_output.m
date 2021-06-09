%% Setup
init_2
plot_result

%% Load Data
load('210520_Versuch');
main_data = data{2};
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


% %% Clean Data (not needed, data is already clean)
% [volt_clean, place] = clean(voltage,8);
% cur_clean = timeseries(current.Data(place), current.Time(place));


%% Find Variance of zeta Z_k
[m,ind] = min(abs(voltage.Time-100));
static_volt = voltage.Data(1:ind);
Z_k = var(static_volt);

for soc_0 = 0:4/30:0.7
    %% 2nd Version of EKF
    % x1 = 1/C; x2 = 1/R; x3 = 1/R0; 
    % x4 = SOC; x5 = U (Voltage over Capacitor)
    % Initial state
    % x1_0 = 1/790.27; x2_0 = 1/0.0126; x3_0 = 1/0.0086;
    % x4_0 = 50%; x5_0 = 0;
    load('converged_param.mat');
    a1 = -4.0568; a2 = 9.7508; a3 = 19.6545;


    uoc = voltage.Data(1)-current.Data(1)*param_conv(3);
    soc = [(-a2+sqrt(a2^2-4*a1*(a3-uoc))),(-a2+sqrt(a2^2-4*a1*(a3-uoc)))]/(2*a1);
    if soc(1)<=1 & soc(1)>=0
        soc = soc(1);
    elseif soc(2)<=1 & soc(2)>=0
        soc = soc(2);
    else
        soc = 0.5;
    end
    Delta = 8e-3; 
    x0 = [1/param_conv(1); 1/param_conv(2); 1/param_conv(3); soc_0; 0];
    Q_e = 2.9*3600;
    P0 = eye(5); P0(1,1)=(0.03*x0(1))^2;P0(2,2)=(0.03*x0(2))^2;P0(3,3)=(0.03*x0(3))^2;
    P0(4,4)=0.5^2;P0(5,5)=(voltage.Data(1)*0.05)^2;

    W_k = 0.0*eye(5);

    %% Simulate
    filename = 'akku_schaetzung';
    systemModel = load_system(filename);
    modelConfigSet = getActiveConfigSet(systemModel);
    modelConfig = modelConfigSet.copy;
    simResult = sim(filename, modelConfig);

    t = simResult.tout';
    y_sim = simResult.simout_check(:,1)';
    u = simResult.simout(:,6);
    y = simResult.simout(:,7)';
    titl = ['Model Validation: Dataset 2, SOC_0: ', num2str(soc_0)];

    figure('Name', titl);
    subplot(3,1,1);
    plot(t, y, 'b', 'LineWidth',1.4);
    hold on
    plot(t, y_sim, 'g--', 'LineWidth',0.8);
    legend({'Measured Output y_{mess}','Simulated Output y_{sim}'},...
        'Location','northwest');
    grid on; hold off
    ylabel('Voltage [V]');
    subplot(3,1,2);
    plot(t, y-y_sim, 'r');
    ylabel('Voltage [V]');
    legend({'Error between y_{mess} and y_{sim}'},...
        'Location', 'northwest');
    grid on;
    subplot(3,1,3);
    plot(t, u, 'r');
    ylabel('Current [A]'); xlabel('Time [s]');
    legend({'Current u'});
    grid on;
    
    %pic_file = ['validation_soc_', num2str(soc_0),'.pdf'];
    %exportgraphics(gcf,pic_file,'ContentType','Vector');
end


% identify SoC-Uoc Curve using system identification toolbox
%clc;%clear;
%% Setup
if 1
    load('src/210609_data');
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
    time_est = 300;
    time_val = time_mes - time_est;

    [~,ind_val] = min(abs(voltage.Time-time_est));

    volt_est = timeseries(voltage.Data(1:ind_val),voltage.Time(1:ind_val));
    volt_val = timeseries(voltage.Data(ind_val+1:end),...
        voltage.Time(ind_val+1:end)-voltage.Time(ind_val+1));

    cur_est = timeseries(current.Data(1:ind_val),current.Time(1:ind_val));
    cur_val = timeseries(current.Data(ind_val+1:end),...
        current.Time(ind_val+1:end)-current.Time(ind_val+1));

    est_data = iddata(volt_est.Data, cur_est.Data, [], 'SamplingInstants',...
        cur_est.Time);
    val_data = iddata(volt_val.Data, cur_val.Data,[],'SamplingInstants',...
        cur_val.Time);

%% Find Variance of zeta Z_k
    [m,ind] = min(abs(voltage.Time-55));
    static_volt = voltage.Data(1:ind);
    Z_k = var(static_volt);

%% SOC-Graph
    u_oc = 2*[9.75, 10.1, 10.3, 10.55, 10.75, 10.9, 11.2, 11.3, 11.45, ...
        11.55, 11.7, 11.8, 12.05, 12.15, 12.25, 12.35, 12.46, 12.55, ...
        12.60, 12.60, 12.65]';
    soc = 0:0.05:1;
    phi = [soc'.^2 soc' ones(21,1)];
    a = phi'*phi\phi'*u_oc; 
    a_ = a;

    figure();
    plot(soc, u_oc); hold on
    plot(soc, phi*a);
    grid on; hold off

%exportgraphics(gcf, 'graph/soc_kurve.pdf', 'ContentType', 'Vector');
%% EKF
% x1 = 1/C; x2 = 1/R; x3 = 1/R0; 
% x4 = SOC; x5 = U (Voltage over Capacitor)
% Initial state
% x1_0 = 1/790.27; x2_0 = 1/0.0126; x3_0 = 1/0.0086;
% x4_0 = 50%; x5_0 = 0;

    Delta = 8e-3; 
    x0 = [2/790.27; 1/0.0126/2; 1/0.0086/2; 0.5; 0];
    Q_e = 2.9*3600;

    P0 = eye(5); P0(1,1)=(0.03*x0(1))^2;P0(2,2)=(0.03*x0(2))^2;P0(3,3)=(0.03*x0(3))^2;
    P0(4,4)=0.5^2;P0(5,5)=(voltage.Data(1)*0.05)^2;

    W_k = 0.0*eye(5);
end

%% Simulate EKF
if 1
    % Activate EKF
    filename = 'akku_schaetzung';
    systemModel = load_system(filename);
    modelConfigSet = getActiveConfigSet(systemModel);
    modelConfig = modelConfigSet. copy;
    simResult = sim(filename, modelConfig);

    x = simResult.simout(1:6,:);
    x_conv = x(:,end);

    % Validate using level 2 s-function
    filename = 'validierung';
    systemModel = load_system(filename);
    modelConfigSet = getActiveConfigSet(systemModel);
    modelConfig = modelConfigSet.copy;
    simResult = sim(filename, modelConfig);
    t_val = simResult.tout;
    y_hat = simResult.simout(:,1);
end
%% simulate with idnlgrey
order = [1, 1, 2];
bat_model = idnlgrey('./function/lead_battery_nl',order,[a_;x_conv(1:3);...
    x_conv(3)],x_conv(4:5));
% fix parameter value of SoC-Uoc Curve
if 0
    for i = 1:3
        bat_model.Parameters(i).Fixed = true;
    end
else
    for i = 4:7
        bat_model.Parameters(i).Fixed = true;
    end
end
    
opt = nlgreyestOptions('Display', 'on');
opt.SearchMethod = 'fmincon';
opt.SearchOptions.Algorithm = 'trust-region-reflective';
opt.SearchOptions.MaxIterations = 20;
est_sqp_nor = nlgreyest(est_data, bat_model, opt);

y_sim = sim(est_sqp_nor, val_data);
%fit = repmat(struct('SearchMethod', '', 'Value', []),5,1);
%fit(5).SearchMethod = 'fmincon-trr with fixed soc';
%fit(6).Value = (1-goodnessOfFit(y_sim.y,val_data.y,'NRMSE'))*100;



figure(); hold on
plot(t_val, y_hat);
plot(val_data.SamplingInstants, val_data.y, 'g');
plot(val_data.SamplingInstants, y_sim.y, '--');
grid on; hold off
legend('Estimated-EKF', 'Data', 'Estimated');



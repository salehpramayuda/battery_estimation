%% Load Data
data_set = 1;
load('210520_Versuch');
main_data = data{data_set};
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
t = voltage.Time;
u = current.Data;
y = voltage.Data;

%% Power plot
p = y; 
for i = 1:length(p)
    p(i) = u(i)*y(i);
end

figure();
subplot(2,1,1);plot(t,u,'r');ylabel('Current [A]');grid on;
title(['Power Progression: Dataset ', num2str(data_set)]);
subplot(2,1,2);plot(t,p,'g');ylabel('Power [W]');xlabel('Time [s]');grid on;
exportgraphics(gcf,'power_verlauf_1.pdf','ContentType','Vector');

%% RMS Power
rms_power = 1/t(end)*sqrt(sum(p.^2));
disp(rms_power);

%% Energy
energy = 0;
e_v = [];
for i = 1:length(p)-1
    energy = p(i)*(t(i+1)-t(i))+ energy;
    e_v = [e_v energy];
end
disp(energy);
figure();plot(t(1:end-1),e_v, 'g');grid on;
ylabel('Harvested Energy [J]');xlabel('Time [s]');
title(['Energy Progression: Dataset ',num2str(data_set)]);
exportgraphics(gcf,'energie_verlauf.pdf','ContentType','Vector');

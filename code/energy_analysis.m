load src/210609_data.mat
t=data{11}{12}.Values.time;
v=data{11}{12}.Values.data;
V=data{11}{13}.Values.data;
I=data{11}{14}.Values.data;
F=data{11}{8}.Values.data;
Ce=195;
P_damping_electric=v.^2*Ce;
%P_damping_electric=v.*F(1:5:end);
P_harvested=(I+0.025).*V;
P=I.*V;

figure(1);clf;
subplot(2,1,1);
plot(t,P_damping_electric,'b'); hold on;
plot(t,P_harvested,'g');
plot(t,-0.025*V,'r');
ylabel('Power in W');
title('Time in s');
legend('P dissipated by Stepper/EHC','P harvested into battery','P consumption MC+Sensors');

figure(2);clf;
plot(P_harvested./P_damping_electric,'b'); hold on;

mean(P_harvested)/mean(P_damping_electric)

E=cumtrapz(t,P)

figure(1);
subplot(2,1,2);
plot(t,E);
grid on;
ylabel('Energy in J');
title('Time in s');







%Design 1
% Fast, large overshoot
Kp = 0.13996;
Kd = 0.0058214;

%Design 2
% Fast, small overshoot
Kp = 0.40385;
Kd = 0.017029;

%Design 3
% Minimize control effort
Kp = 0.064719;
Kd =  0.017472;

 

figure
sim('pdstep2.slx')
plot(simout.Time, simout.Data)
hold on
pddata = importdata("Step_PD_Kp0.064719_Kd0.017472_1.csv");
plot(pddata.data(:,1), pddata.data(:, 10))
xlim([0, 2])

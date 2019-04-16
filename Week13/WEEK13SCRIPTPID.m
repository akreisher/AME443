
%Design 1
%Decent speed/CE, small overshoot
Kp = 0.33267;
Kd = 0.017472;
Ki = 1.7837;

%Design 2
% Minimize ce
Kp = 0.05492;
Kd = 0.017472;
Ki = 0.40662;

%Design 3
% Minimize overshoot
Kp = 0.23487;
Kd = 0.017472;
Ki = 0.95892;


figure
sim('pidstep2.slx')
plot(simout.Time, simout.Data)
hold on
pddata = importdata("Vel FB/Step_PID_Kp0.23487_Kd0.017472_Ki0.95892_1.csv");
ce = pddata.data(:, 5)*5/16383;
plot(pddata.data(:,1), pddata.data(:, 10))

xlim([0, 2])

figure
plot(pddata.data(:,1), ce)
hold on
plot(cesim.Time, cesim.Data/100)
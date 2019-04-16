%Design 1
% Small overshoot ss error
Kp = 0.57273;
Kd = 0.019504;

%Design 2
% Quicker
%  Kp = 0.086231;
%  Kd = 0.006507;



figure
sim('pdstepserial.slx')
plot(simout.Time, simout.Data)
hold on
pddata = importdata("Step_PD_Kp0.57273_Kd0.019504_1.csv");
ce = pddata.data(:, 5)*5/16383;
plot(pddata.data(:,1), pddata.data(:, 10))
xlim([0, 2])

figure()
plot(pddata.data(:,1), ce)

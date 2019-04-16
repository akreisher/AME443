%Design 1
%  minimize overshoot

Kp = 0.1142;
Kd = 0.0081717;
Ki = 0.50756;

%Design 2
%  Quicker settling

Kp = 0.0604;
Kd = 0.0075369;
Ki = 0.97295;


figure
sim('pdstepserial.slx')
plot(simout.Time, simout.Data)
hold on
pddata = importdata("Step_PID_Kp0.1142_Kd0.0081717_Ki0.50756_1.csv");
ce = pddata.data(:, 5)*5/16383;
plot(pddata.data(:,1), pddata.data(:, 10))
xlim([0, 3])

figure()
plot(pddata.data(:,1), ce)
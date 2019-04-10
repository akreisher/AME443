data = importdata("StepCL_r100_1.csv");

t = data.data(:,1);
pos1 = data.data(:, 10);
pos2 = data.data(:, 14);

t = t(t<2);
pos1 = pos1(t<2);
pos2 = pos2(t<2);

ss1 = pos1(end);

omega1 = 17;
omegaz = 30;
omega2 = 43;

zeta1 = 0.048/2;
zeta2 = 0.05/2;
zetaz = 0.163/2;

A = 100;

KHWCL = ss1*((omega1*omega2/omegaz)^2)/A;


numcol = KHWCL*[1, 2*zetaz*omegaz, omegaz^2];
dencol = conv([1, 2*zeta1*omega1, omega1^2], [1, 2*zeta2*omega2, omega2^2]);

sys = tf(num, den);

sim('step.slx')
plot(simout.Time, simout.Data)
hold on
plot(t, pos1)



% PD Control
% Design 1 - pretty quick, small overshoot/ce
% Kd = 0.041702;
% Kp = 1.0428;

%design 2 - no overshoot, minimal ce, kinda slow
Kd = 0.022673;
Kp = 0.10636;

%design 3 - Fast, small overshoot, largeish CE
% Kd = 0.082453;
% Kp =  2.6615;

figure
sim('pdstep.slx')
plot(simout.Time, simout.Data)
hold on
pddata = importdata("Step_PD_Kp0.10636_Kd0.022673_1.csv");
plot(pddata.data(:,1), pddata.data(:, 10))
xlim([0, 2])


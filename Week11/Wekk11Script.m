data = importdata("SineSweep_2DOF_m1484_m2_485_k3k3_r0.2_e10_2.csv");

t = data.data(:,1);
pos1 = data.data(:,10);
pos2 = data.data(:, 14);
ce = data.data(:, 5);

[pks, ind] = findpeaks(ce);

T = diff(t(ind));

omegas = 2*pi./T;

ts = t(ind(2:end));

f = fit(ts, omegas, 'poly1');

a = f.p1;
b = f.p2;

[pks,ind] = findpeaks(pos1);
dB1 = 20*log10(abs(pks./0.2));
omegas = a*t(ind)+b;

[pks,ind] = findpeaks(pos2);
dB2 = 20*log10(abs(pks./0.2));

figure
plot(omegas, dB1, 'Color', [0 0.4470 0.7410] )
hold on
plot(omegas, dB2(1:length(omegas)), 'Color', [0.8500 0.3250 0.0980])

files = (dir("Sine"));
files = files(3:end);

singleomegas = zeros(1, length(files));
singlegains1 = zeros(1, length(files));
singlegains2 = zeros(1, length(files));
for i = 1:length(files)
    data = importdata("Sine/"+files(i).name);
    j = strfind(files(i).name,'r');
    x = str2num(files(i).name(j+1:j+3));
    t = data.data(:,1);
    pos1 = data.data(:,10);
    pos2 = data.data(:, 14);
    ce = data.data(:, 5);
    
    [pks, ind] = findpeaks(ce);
    T = diff(t(ind));

    singleomegas(i) = mean(2*pi./T);
    
    [pks, ind] = findpeaks(pos1);
    singlegains1(i) = 20*log10(pks(end)/x);
    
    [pks, ind] = findpeaks(pos2);
    singlegains2(i) = 20*log10(pks(end)/x);
    
end
[singleomegas, ind] = sort(singleomegas);
plot(singleomegas, singlegains1(ind),'o', 'Color', [0 0.4470 0.7410]);
hold on
plot(singleomegas, singlegains2(ind),'o', 'Color', [0.8500 0.3250 0.0980]);
set(gca, 'xscale', 'log')



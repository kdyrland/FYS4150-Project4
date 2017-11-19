%% Project 4

% part b
load('proj4.mat')

numE(1:length(mcycl)) = -7.983928;  % analytical value
nummM(1:length(mcycl)) = 3.994643;

% plot
figure
semilogx(mcycl,energy/2^2,'Displayname','Numerical')
hold on
semilogx(mcycl,numE/2^2,'Displayname','Analytical')
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<E>/L^2')
lgd = legend('show');
lgd.FontSize = 14;

% figure
% plot(mcycl,abs(mag))
% xlabel('Monte Carlo cycles','Fontsize',14)
% ylabel('<M>')

figure
semilogx(mcycl,meanmag/2^2,'Displayname','Numerical')
hold on
semilogx(mcycl,nummM/2^2,'Displayname','Analytical')
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<|M|>/L^2')
lgd = legend('show');
lgd.FontSize = 14;

%% part c

load('proj4.mat')

% accepted configs
mc1 = mcycac(1:19);
mc2 = mcycac(20:end);
ac1 = aconfig(1:19);
ac2 = aconfig(20:end);
ac1n = aconfig(1:19)./max(aconfig(1:19));
ac2n = aconfig(20:end)./max(aconfig(20:end));

figure
loglog(mc1,ac1,'Displayname','T = 1 J/k_B')
hold on
loglog(mc2,ac2,'Displayname','T = 2.4 J/k_B')
ylabel('Accepted configurations','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
lgd = legend('show','Location','southeast');
lgd.FontSize = 14;

% figure
% hold on
% plot(mc1,ac1n,'Displayname','T = 1 J/k_B')
% plot(mc2,ac2n,'Displayname','T = 2.4 J/k_B')
% ylabel('Accepted configurations (normalized)','Fontsize',14)
% xlabel('Monte Carlo cycles','Fontsize',14)

% expectation values
% T = 1.0 - ordered
c_mc = partc1ord(:,1);  % mc cycles
c_e1o = partc1ord(:,2); % energy
c_m1o = partc1ord(:,3); % mag

figure
semilogx(c_mc,c_e1o./20^2)
title('T = 1.0 J/k_B','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<E>/L^2','Fontsize',16)

figure
semilogx(c_mc,c_m1o./20^2)
title('T = 1.0 J/k_B','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<|M|>/L^2','Fontsize',16)


% T = 1.0 - random
c_e1r = partc1rnd(:,2);     % energy
c_m1r = partc1rnd(:,3);     % mag

figure
semilogx(c_mc,c_e1r/20^2)
title('T = 1.0 J/k_B','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<E>/L^2','Fontsize',16)

figure
semilogx(c_mc,abs(c_m1r)/20^2)
title('T = 1.0 J/k_B','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<|M|>/L^2','Fontsize',16)


% T = 2.4 - ordered
c_e24o = partc24ord(:,2);   % energy
c_m24o = partc24ord(:,3);   % mag

figure
semilogx(c_mc,c_e24o./20^2)
title('T = 2.4 J/k_B','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<E>/L^2','Fontsize',16)

figure
semilogx(c_mc,c_m24o/20^2)
title('T = 2.4 J/k_B','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<M>/L^2','Fontsize',16)


% T = 2.4 - random
c_e24r = partc24rnd(:,2);   % energy
c_m24r = partc24rnd(:,3);   % mag

figure
semilogx(c_mc,c_e24r/20^2)
title('T = 2.4 J/k_B','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<E>/L^2','Fontsize',16)

figure
semilogx(c_mc,abs(c_m24r)/20^2)
title('T = 2.4 J/k_B','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<M>/L^2','Fontsize',16)


%% part d

load('proj4.mat')

% plot histogram
bins = 40;  % histogram bins

figure
hold on
% plot histogram of last value in each random walk
histogram(nrj1,bins,'Normalization','probability','Displayname','Numerical')
title('Probability Distribution for T = 1.0 J/k_B','Fontsize',14)
xlabel('<E>','Fontsize',16)
ylabel('Normalized Probabilty','Fontsize',14)

figure
hold on
% plot histogram of last value in each random walk
histogram(nrj24,bins,'Normalization','probability','Displayname','Numerical')
title('Probability Distribution for T = 2.4 J/k_B','Fontsize',14)
xlabel('<E>/L^2','Fontsize',16)
ylabel('Normalized Probabilty','Fontsize',14)



%% part e

%energy
figure
hold on
plot(etemp,nrj20,'Displayname','L = 20')
plot(etemp,nrj40,'Displayname','L = 40')
plot(etemp,nrj60,'Displayname','L = 60')
plot(etemp,nrj80,'Displayname','L = 80')
plot(etemp,nrj100,'Displayname','L = 100')
xlabel('Temperature [J/k_B]','Fontsize',14)
ylabel('<E>/L^2','Fontsize',16)

lgd = legend('show','Location','southeast');
lgd.FontSize = 14;

%% specific heat
figure
hold on
plot(etemp,cv20,'Displayname','L = 20')
plot(etemp,cv40,'Displayname','L = 40')
plot(etemp,cv60,'Displayname','L = 60')
plot(etemp,cv80,'Displayname','L = 80')
plot(etemp,cv100,'Displayname','L = 100')
xlabel('Temperature [J/k_B]','Fontsize',14)
ylabel('C_v','Fontsize',16)

lgd = legend('show');
lgd.FontSize = 14;



%% mean magnetization
figure
hold on
plot(etemp,meanmag20,'Displayname','L = 20')
plot(etemp,meanmag40,'Displayname','L = 40')
plot(etemp,meanmag60,'Displayname','L = 60')
plot(etemp,meanmag80,'Displayname','L = 80')
plot(etemp,meanmag100,'Displayname','L = 100')
xlabel('Temperature [J/k_B]','Fontsize',14)
ylabel('<|M|>','Fontsize',16)

lgd = legend('show');
lgd.FontSize = 14;

%% susceptibility
figure
hold on
plot(etemp,xsi20,'Displayname','L = 20')
plot(etemp,xsi40,'Displayname','L = 40')
plot(etemp,xsi60,'Displayname','L = 60')
plot(etemp,xsi80,'Displayname','L = 80')
plot(etemp,xsi100,'Displayname','L = 100')
xlabel('Temperature [J/k_B]','Fontsize',14)
ylabel('\chi','Fontsize',20)

lgd = legend('show');
lgd.FontSize = 14;






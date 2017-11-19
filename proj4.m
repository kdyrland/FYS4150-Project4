%% Project 4

% part b
load('proj4.mat')

numE(1:length(mcycl)) = -7.983928;  % analytical value
nummM(1:length(mcycl)) = 3.994643;

%% plot
figure
plot(mcycl,energy/2^2,'Displayname','Numerical')
hold on
plot(mcycl,numE/2^2,'Displayname','Analytical')
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<E>/L^2','Fontsize',16)
xlim([-1e5 1e6])
lgd = legend('show');
lgd.FontSize = 14;

% figure
% plot(mcycl,abs(mag))
% xlabel('Monte Carlo cycles','Fontsize',14)
% ylabel('<M>')

figure
plot(mcycl,meanmag/2^2,'Displayname','Numerical')
hold on
plot(mcycl,nummM/2^2,'Displayname','Analytical')
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<|M|>/L^2','Fontsize',16)
xlim([-1e5 1e6])
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

%% expectation values
% T = 1.0 - ordered
c_mc = partc1ord(:,1);  % mc cycles
c_e1o = partc1ord(:,2); % energy
c_m1o = partc1ord(:,3); % mag

figure
plot(c_mc,c_e1o./20^2)
title('T = 1.0 J/k_B','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<E>/L^2','Fontsize',16)
xlim([-1e5 1e6])

figure
plot(c_mc,c_m1o./20^2)
title('T = 1.0 J/k_B','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<|M|>/L^2','Fontsize',16)
xlim([-1e5 1e6])


% T = 1.0 - random
c_e1r = partc1rnd(:,2);     % energy
c_m1r = partc1rnd(:,3);     % mag

figure
plot(c_mc,c_e1r/20^2)
title('T = 1.0 J/k_B','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<E>/L^2','Fontsize',16)
xlim([-1e5 1e6])
ylim([1.3 2.1])

figure
plot(c_mc,abs(c_m1r)/20^2)
title('T = 1.0 J/k_B','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<|M|>/L^2','Fontsize',16)
xlim([-1e5 1e6])
ylim([0.2 1.1])


% T = 2.4 - ordered
c_e24o = partc24ord(:,2);   % energy
c_m24o = partc24ord(:,3);   % mag

figure
plot(c_mc,c_e24o./20^2)
title('T = 2.4 J/k_B','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<E>/L^2','Fontsize',16)
xlim([-1e5 1e6])

figure
plot(c_mc,c_m24o/20^2)
title('T = 2.4 J/k_B','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<M>/L^2','Fontsize',16)
xlim([-1e5 1e6])


% T = 2.4 - random
c_e24r = partc24rnd(:,2);   % energy
c_m24r = partc24rnd(:,3);   % mag

figure
plot(c_mc,c_e24r/20^2)
title('T = 2.4 J/k_B','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<E>/L^2','Fontsize',16)
xlim([-1e5 1e6])

figure
plot(c_mc,abs(c_m24r)/20^2)
title('T = 2.4 J/k_B','Fontsize',14)
xlabel('Monte Carlo cycles','Fontsize',14)
ylabel('<M>/L^2','Fontsize',16)
xlim([-1e5 1e6])


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


%% Find Tc
% using Cv
maxCv20 = find(cv20 == max(cv20));
maxCv40 = find(cv40 == max(cv40));
maxCv60 = find(cv60 == max(cv60));
maxCv80 = find(cv80 == max(cv80));
maxCv100 = find(cv100 == max(cv100));

Tc_cv20 = etemp(maxCv20);
Tc_cv40 = etemp(maxCv40);
Tc_cv60 = etemp(maxCv60);
Tc_cv80 = etemp(maxCv80);
Tc_cv100 = etemp(maxCv100);
disp(['T_C(L=20) = ', num2str(Tc_cv20)])
disp(['T_C(L=40) = ', num2str(Tc_cv40)])
disp(['T_C(L=60) = ', num2str(Tc_cv60)])
disp(['T_C(L=80) = ', num2str(Tc_cv80)])
disp(['T_C(L=100) = ', num2str(Tc_cv100)])
disp(newline)

% using chi
maxchi20 = find(xsi20 == max(xsi20));
maxchi40 = find(xsi40 == max(xsi40));
maxchi60 = find(xsi60 == max(xsi60));
maxchi80 = find(xsi80 == max(xsi80));
maxchi100 = find(xsi100 == max(xsi100));

Tc_chi20 = etemp(maxchi20);
Tc_chi40 = etemp(maxchi40);
Tc_chi60 = etemp(maxchi60);
Tc_chi80 = etemp(maxchi80);
Tc_chi100 = etemp(maxchi100);
disp(['T_C(L=20) = ', num2str(Tc_chi20)])
disp(['T_C(L=40) = ', num2str(Tc_chi40)])
disp(['T_C(L=60) = ', num2str(Tc_chi60)])
disp(['T_C(L=80) = ', num2str(Tc_chi80)])
disp(['T_C(L=100) = ', num2str(Tc_chi100)])
disp(newline)


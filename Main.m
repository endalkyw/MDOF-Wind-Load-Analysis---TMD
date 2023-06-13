%% The latest version of the project without UQ
clear; 
clc;
run("4. Matlab Scripts\func_coll(""run_"")");

%% Displacement Plots 
di = readmatrix("2. Outputs\M2 RelativeDisplacement_TH.csv");
tm = readmatrix("2. Outputs\M1 TimeVector_TH.csv");

%
hold on
plot(tm,di(18,:),'r--');
grid on
xlabel("time(s)");
ylabel("tip displacement(m)");

%%
wl = readmatrix("3. Wind TH\WindLoad.csv");
ws = readmatrix("3. Wind TH\WindSpeed.csv");
% wl = wl(:,1:6000);
% writematrix(wl,"3. Wind TH\WindLoad.csv")
hold on
plot(wl(18,:))

%%
d1 = readmatrix("2. Outputs\sim_1\M2 RelativeDisplacement_TH.csv");
d2 = readmatrix("2. Outputs\sim_2\M2 RelativeDisplacement_TH.csv");
d3 = readmatrix("2. Outputs\sim_3\M2 RelativeDisplacement_TH.csv");
%
hold on
plot(d1(18,:),'b-');
plot(d2(18,:),'r-');
plot(d3(18,:),'k-');
hold off




















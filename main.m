clear all
close all
clc

%% Model parameters
% *1: da "Hybrid Electric Energy Management Strategies, Onori, Serrao, Rizzoni"
parameters.rhoa     =   1.22;   % *1 [kg/m^3]
parameters.Af       =   2.33;   % *1 [m^2]
parameters.Cd       =   0.26;   % *1 
parameters.g        =   9.81;   % [N*m/s^2]
parameters.Crr      =   0.024;  % *1 
parameters.m_v      =   1370;   % *1 [kg]
parameters.R        =   0.32;   % *1 [m]
parameters.lambda   =   43.4*1e6;
parameters.eta_em   =   0.8;
parameters.Voc      =   240;    % *1 [V]
parameters.R0       =   0.13;   % *1 [Ohm]
parameters.Qnom     =   23400;  % *1 [C]
parameters.eta_coul =   0.95;

parameters.reg_break_limit = 10000; % regenerative breaking max recharge power [kW]

dim_decision_var = 4; % time interval in seconds, equivalent to the dimension of the decision variable
interval_size = 1800/dim_decision_var; % will be length of u_compressed
parameters.interval_size = interval_size;

% lower and higher constraints on soc, 0.1 < soc_vec(t) < 0.9 for every t
parameters.soc_lowerConstraint = 0.1;
parameters.soc_higherConstraint = 0.9;
parameters.initial_soc = 0.8; % initial state of charge of the battery

% read wltp excel file
wltp_cycle = xlsread('wltp_cycle.xlsx');
parameters.time_vec = wltp_cycle(:,1);
parameters.v_vec = wltp_cycle(:,2)/3.6;
parameters.a_vec = wltp_cycle(:,3);

%% Plot speed and acceleration profile
% figure(1),
% subplot(2,1,1),plot(parameters.time_vec, parameters.v_vec*3.6),grid on, hold on,xlabel('time [s]'),ylabel('v_ref [km/h]')
% subplot(2,1,2),plot(parameters.time_vec, parameters.a_vec),grid on, hold on,xlabel('time [s]'),ylabel('a_ref [m/s^2]')

%% Retrive efficiency map
load lasso_param.mat
parameters.lasso_param=lasso_param;

%% Simulation parameters
Ts=1; % sampling time (s)
Tend=1800; % end time (s)
tvec=parameters.time_vec;
N=length(tvec);
parameters.N = N;

%% Linear equality constraint parameters
% there are no linear equality constraints
A = [];
b = [];

%% Linear inequality constraint parameters
% constraint u<=1 written as -u >= -1, C*u >= d
C = -eye(interval_size);
d = -ones(interval_size,1);

%% Optimization: Gauss Newton
% Initialize solver options
myoptions               =   myoptimset;
myoptions.Hessmethod  	=	'GN';
myoptions.gradmethod  	=	'FD';
myoptions.graddx        =	2^-26;
myoptions.tolgrad    	=	1e-10;
myoptions.ls_tkmax      =	1;          
myoptions.ls_beta       =	0.1;
myoptions.ls_c          =	0.1;
myoptions.ls_nitermax   =	10;

myoptions.nitermax      =	1;

myoptions.xsequence     =	'on';
myoptions.GN_funF       =   @(x)big_fun(x,parameters);
myoptions.GN_sigma      =   1e-3;

% set initial state x0
x0 = zeros(interval_size,1);
% declare the number of nonlinear constraints
p = 1; % 1 equality constraint defined as "g = soc_vec(N,1) - soc_vec(1,1);" in "big_fun.m" and "fun.m"
q = 2*N; % N inequality constraints "soc_vec(ind,1) > 0.1 for every Ts" and N inequality constraints "soc_vec(ind,1) < 0.9 for every Ts"

%% Run solver
tic
[xstar,fxstar,niter,exitflag,xsequence] = myfmincon(@(x)fun(x,parameters),x0,A,b,C,d,p,q,myoptions);
toc
%% xstar_decompressed computation:
xstar_decompressed=[];
for i=1:length(xstar)
    for j=1:round(N/interval_size)
        xstar_decompressed = [xstar_decompressed ; xstar(i)];
    end
end
xstar_decompressed = [xstar_decompressed; xstar_decompressed(end)];

%% Simulation (only ICE)
m_f_ICE= zeros(N,1); % mass of the fuel at every second
soc_vec_ICE = zeros(N,1); % state of charge at every second

% initialization of the soc
soc_vec_ICE(1,1) = parameters.initial_soc;

% decision variable when we only use ICE is u(t)=0 for every t
u_vec = zeros(N,1);

% simulation
for ind=2:N
    [m_f_dot,soc_dot] = fuel_consumption(parameters, u_vec(ind,1),ind);
    m_f_ICE(ind,1) = m_f_ICE(ind-1,1)+Ts*m_f_dot; 
    soc_vec_ICE(ind,1) = soc_vec_ICE(ind-1,1); % soc never changes
end

%% Simulation ICE only + EM from regenerative braking
m_f_ICEreg = zeros(N,1);
soc_vec_ICEreg = zeros(N,1);

% initialization of the soc
soc_vec_ICEreg(1,1) = parameters.initial_soc;

% second strategy -> use EM if soc > initial_soc
u_vec = zeros(N,1);

for ind=2:N
    [m_f_dot,soc_dot] = fuel_consumption(parameters, u_vec(ind,1),ind);
    
    m_f_ICEreg(ind,1) = m_f_ICEreg(ind-1,1)+Ts*m_f_dot; 

    soc_vec_ICEreg(ind,1) = soc_vec_ICEreg(ind-1,1)+Ts*soc_dot;
    
    if soc_vec_ICEreg(ind,1) > parameters.initial_soc
        u_vec(ind+1,1) = 1;
    end
end
%% Simulation (optimum control)
m_f_opt = zeros(N,1); % mass of the fuel
soc_vec_opt = zeros(N,1); % state of charge

% initialization of the soc -> start from a safe segion SOC=0.9
soc_vec_opt(1,1) = parameters.initial_soc;

% compute xstar_decompressed
xstar_decompressed=[];
for i=1:length(xstar)
    for j=1:round(N/interval_size)
        xstar_decompressed = [xstar_decompressed ; xstar(i)];
    end
end
xstar_decompressed = [xstar_decompressed; xstar_decompressed(end)];

% simulation
u_vec = xstar_decompressed;
for ind=2:N
    [m_f_dot,soc_dot] = fuel_consumption(parameters, u_vec(ind,1),ind);
    m_f_opt(ind,1) = m_f_opt(ind-1,1)+Ts*m_f_dot;
    soc_vec_opt(ind,1) = soc_vec_opt(ind-1,1)+Ts*soc_dot;
end

%% Plot the results
close all

figure(1),
plot(tvec, m_f_ICE,'LineWidth',1.5),grid on, hold on,xlabel('time [s]'),ylabel('fuel consumption [kg]')
% plot(tvec, m_f_ICEreg,'LineWidth',1.5,'color', 'r'),grid on, hold on,xlabel('time [s]'),ylabel('fuel consumption [kg]')
plot(tvec, m_f_opt,'LineWidth',1.5,'color', 'g'),grid on, hold on,xlabel('time [s]'),ylabel('fuel consumption [kg]')
legend('ICE only', 'ICE+RB', 'Optimum Ctrl')

figure(2),
plot(tvec, soc_vec_ICE,'LineWidth',1.5),grid on, hold on, xlabel('time [s]'),ylabel('soc')
% plot(tvec, soc_vec_ICEreg,'LineWidth',1.5,'color', 'r'),grid on, hold on, xlabel('time [s]'),ylabel('soc')
plot(tvec, soc_vec_opt,'LineWidth',1.5,'color', 'g'),grid on, hold on, xlabel('time [s]'),ylabel('soc')
legend('ICE only', 'ICE+RB','Optimum Ctrl')
axis([-inf inf 0 1])

figure(3)
plot(tvec, xstar_decompressed,'LineWidth',1.5),grid on, hold on, xlabel('time [s]'),ylabel('control variable')
legend('Optimum Ctrl')
axis([-inf inf -2 1.2])

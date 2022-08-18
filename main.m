clear all
close all
clc

%% Model parameters
% *1: da "Hybrid Electric Energy Management Strategies, Onori, Serrao, Rizzoni"
parameters.rhoa     =   1.22;   % *1
parameters.Af       =   2.33;   % *1
parameters.Cd       =   0.26;   % *1
parameters.g        =   9.81;   
parameters.Crr      =   0.024;  % *1
parameters.m_v      =   1370;   % *1
parameters.R        =   0.32;   % *1
parameters.lambda   =   43.4*1e6;
parameters.eta_em   =   0.8;
parameters.Voc      =   240;    % *1 should be variable wrt state of charge, see pag.92
parameters.R0       =   0.13;   % *1 should be variable wrt state of charge, see pag.92
parameters.Qnom     =   23400;  % *1
parameters.eta_coul =   0.95;

dim_decision_var = 4;
parameters.interval_size = 1800/dim_decision_var; % (dimension of the decision variable) will be length of u_compressed

parameters.initial_soc = 0.8;

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
Ts=1;            % sampling time (s)
Tend=1800;
tvec=parameters.time_vec;
N=length(tvec);
parameters.N = N;
interval_size = parameters.interval_size; % (dimension of the decision variable) will be length of u_compressed

%% Linear equality constraint parameters
A = [];
b = [];

%% Linear inequality constraint parameters
% -u >= -1 therefore u<=1
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

myoptions.nitermax      =	5;

myoptions.xsequence     =	'on';
myoptions.GN_funF       =   @(x)big_fun(x,parameters);
myoptions.GN_sigma      =   1e-3;

% set x0
x0 = zeros(interval_size,1);
% number of nonlinear constraints
p = 1; % equaliti constrain defined as "g = soc_vec(N,1) - soc_vec(1,1);" in "big_fun.m" and "fun.m"
q = 2*N; % N inequality constraint "soc_vec(ind,1) > 0.1 for every Ts" and N inequality constraint "soc_vec(ind,1) < 0.9 for every Ts"
% Run solver
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
m_f1 = zeros(N,1);        % this will be the integral in tvec of m_f_dot
m_f_dot_vec = zeros(N,1);
soc_dot_vec = zeros(N,1);
soc_vec1 = zeros(N,1);

% initialization of the soc -> start from a safe segion SOC=0.9
soc_vec1(1,1) = parameters.initial_soc;

% first strategy -> only ICE
u_vec = zeros(N,1);


for ind=2:N
    if soc_vec1(ind-1,1) > 0.9
        u_vec(ind,1) = 1;
    end
    [m_f_dot,soc_dot] = fuel_consumption(parameters, u_vec(ind,1),ind);
    m_f_dot_vec(ind,1) = m_f_dot; 
    m_f1(ind,1) = m_f1(ind-1,1)+Ts*m_f_dot; 

    soc_dot_vec(ind,1) = soc_dot;
    soc_vec1(ind,1) = soc_vec1(ind-1,1)+Ts*soc_dot;
end

%% Simulation ICE only + EM from regenerative braking
% m_f2 = zeros(N,1);        % this will be the integral in tvec of m_f_dot
% m_f_dot_vec = zeros(N,1);
% soc_dot_vec = zeros(N,1);
% soc_vec2 = zeros(N,1);
% 
% % initialization of the soc -> start from a safe segion SOC=0.9
% soc_vec2(1,1) = parameters.initial_soc;
% 
% % second strategy -> only EM until we have battery, then switch to ICE  (uncomment and comment prevous line)
% u_vec = zeros(N,1);
% 
% 
% for ind=2:N
%     [m_f_dot,soc_dot] = fuel_consumption(parameters, u_vec(ind,1),ind);
%     
%     m_f_dot_vec(ind,1) = m_f_dot; 
%     m_f2(ind,1) = m_f2(ind-1,1)+Ts*m_f_dot; 
% 
%     soc_dot_vec(ind,1) = soc_dot;
%     soc_vec2(ind,1) = soc_vec2(ind-1,1)+Ts*soc_dot;
%     
%     if soc_vec2(ind,1) <= 0.1
%         u_vec(ind+1,1) = 0;
%     end
% end
%% Simulation (optimum control)
m_f3 = zeros(N,1); % this will be the integral in tvec of m_f_dot
m_f_dot_vec = zeros(N,1);
soc_dot_vec = zeros(N,1);
eta_eng_vec = zeros(N,1);
Pf_vec = zeros(N,1);
Pem_vec = zeros(N,1);
Peng_vec = zeros(N,1);

soc_vec3 = zeros(N,1);

% initialization of the soc -> start from a safe segion SOC=0.9
soc_vec3(1,1) = parameters.initial_soc;

xstar_decompressed=[];
for i=1:length(xstar)
    for j=1:round(N/interval_size)
        xstar_decompressed = [xstar_decompressed ; xstar(i)];
    end
end
xstar_decompressed = [xstar_decompressed; xstar_decompressed(end)];

u_vec = xstar_decompressed;
for ind=2:N
    [m_f_dot,soc_dot] = fuel_consumption(parameters, u_vec(ind,1),ind);
    
    m_f_dot_vec(ind,1) = m_f_dot; 
    m_f3(ind,1) = m_f3(ind-1,1)+Ts*m_f_dot; 

    soc_dot_vec(ind,1) = soc_dot;
    soc_vec3(ind,1) = soc_vec3(ind-1,1)+Ts*soc_dot;
end

%% Plot the results
close all

figure(1),
plot(tvec, m_f1,'LineWidth',1.5),grid on, hold on,xlabel('time [s]'),ylabel('fuel consumption [kg]')
plot(tvec, m_f3,'LineWidth',1.5,'color', 'g'),grid on, hold on,xlabel('time [s]'),ylabel('fuel consumption [kg]')
legend('ICE only', 'Optimum Ctrl')

figure(2),
plot(tvec, soc_vec1,'LineWidth',1.5),grid on, hold on, xlabel('time [s]'),ylabel('soc')
plot(tvec, soc_vec3,'LineWidth',1.5,'color', 'g'),grid on, hold on, xlabel('time [s]'),ylabel('soc')
legend('ICE only', 'Optimum Ctrl')
axis([-inf inf 0 1])

figure(3)
plot(tvec, xstar_decompressed,'LineWidth',1.5),grid on, hold on, xlabel('time [s]'),ylabel('control variable')
legend('Optimum Ctrl')
axis([-inf inf -1.2 1.2])





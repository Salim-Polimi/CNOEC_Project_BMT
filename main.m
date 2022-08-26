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
parameters.initial_soc = 0.8;

%% read wltp excel file
wltp_cycle = xlsread('wltp_cycle2.xlsx');
parameters.time_vec = wltp_cycle(:,1);
parameters.v_vec = wltp_cycle(:,2)/3.6;
parameters.a_vec = wltp_cycle(:,3);

%% Plot speed and acceleration profile from wltp
figure(1), set(gcf, 'Color', 'White'),
sp(1) = subplot(2,1,1);
hold on, grid on, box on
plot(parameters.time_vec, parameters.v_vec*3.6, 'LineWidth',1.5,'Color','b')
xlabel('Time [s]'), ylabel('Speed [km/h]')
sp(2) = subplot(2,1,2);
hold on, grid on, box on
plot(parameters.time_vec, parameters.a_vec, 'LineWidth',1.5,'Color','b')
xlabel('Time [s]'), ylabel('Acceleration [m/s^2]')
linkaxes(sp,'x'), clear sp

%% Retrive efficiency map
load lasso_param.mat
parameters.lasso_param=lasso_param;

%% Simulation parameters
Ts=1;            % sampling time (s)
tvec=parameters.time_vec;
N=length(tvec);
parameters.N = N;
parameters.dim_decision_var = 1;
parameters.interval_size = N/parameters.dim_decision_var; % (dimension of the decision variable) will be length of u_compressed
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

%% Code generation
% codegen big_fun -args {coder.typeof(x0),coder.typeof(parameters)}
% codegen fun -args {coder.typeof(x0),coder.typeof(parameters)}
% codegen fuel_consumption -args {coder.typeof(parameters),coder.typeof(tvec),coder.typeof(p)}

%% Run optimization
tic
[xstar,fxstar,niter,exitflag,xsequence] = myfmincon(@(x)fun(x,parameters),x0,A,b,C,d,p,q,myoptions);
toc

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
    
    [m_f_dot,soc_dot] = fuel_consumption(parameters, u_vec(ind,1),ind);
    m_f_dot_vec(ind,1) = m_f_dot; 
    m_f1(ind,1) = m_f1(ind-1,1)+Ts*m_f_dot; 
    
    if (soc_dot > 0)
       soc_dot=0; 
    end

    soc_dot_vec(ind,1) = soc_dot;
    soc_vec1(ind,1) = soc_vec1(ind-1,1)+Ts*soc_dot;
end

%% Simulation Heuristic
m_f2 = zeros(N,1);        % this will be the integral in tvec of m_f_dot
m_f_dot_vec = zeros(N,1);
soc_dot_vec = zeros(N,1);
soc_vec2 = zeros(N,1);

% initialization of the soc -> start from a safe segion SOC=0.9
soc_vec2(1,1) = parameters.initial_soc;

% heuristic
u_vec = zeros(N,1);

for ind=2:N
    if parameters.v_vec(ind)*3.6 < 30
        u_vec(ind,1) = 1;
    end
    if parameters.v_vec(ind)*3.6 > 80
        u_vec(ind,1) = -0.2;
    end
    
    [m_f_dot,soc_dot] = fuel_consumption(parameters, u_vec(ind,1),ind);
    
    m_f_dot_vec(ind,1) = m_f_dot; 
    m_f2(ind,1) = m_f2(ind-1,1)+Ts*m_f_dot; 

    soc_dot_vec(ind,1) = soc_dot;
    soc_vec2(ind,1) = soc_vec2(ind-1,1)+Ts*soc_dot;
    
end
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

u_vec = xstar_decompressed;
for ind=2:N
    [m_f_dot,soc_dot] = fuel_consumption(parameters, u_vec(ind,1),ind);
    
    m_f_dot_vec(ind,1) = m_f_dot; 
    m_f3(ind,1) = m_f3(ind-1,1)+Ts*m_f_dot; 

    soc_dot_vec(ind,1) = soc_dot;
    soc_vec3(ind,1) = soc_vec3(ind-1,1)+Ts*soc_dot;
end

%% Plot the results
figure(2), set(gcf, 'Color', 'White'),
grid on, hold on,xlabel('time [s]'),ylabel('fuel consumption [kg]')
plot(tvec, m_f1,'LineWidth',1.5,'color', 'b', 'DisplayName', 'ICE Only')
%plot(tvec, m_f2,'LineWidth',1.5,'color', 'r', 'DisplayName', 'Heuristic')
plot(tvec, m_f3,'LineWidth',1.5,'color', 'g', 'DisplayName', 'Optimum Ctrl')
legend show

figure(3), set(gcf, 'Color', 'White'),
grid on, hold on, xlabel('time [s]'),ylabel('soc [-]'), axis([-inf inf 0 1])
plot(tvec, soc_vec1,'LineWidth',1.5,'color', 'b', 'DisplayName', 'ICE Only')
%plot(tvec, soc_vec2,'LineWidth',1.5,'color', 'r', 'DisplayName', 'Heuristic')
plot(tvec, soc_vec3,'LineWidth',1.5,'color', 'g', 'DisplayName', 'Optimum Ctrl')
legend show

figure(4), set(gcf, 'Color', 'White'),
sp(1) = subplot(3,1,1);
hold on, grid on, box on
plot(parameters.time_vec, parameters.v_vec*3.6, 'LineWidth',1.5,'Color','b')
xlabel('Time [s]'), ylabel('Speed in the driving cycle [km/h]')
sp(2) = subplot(3,1,2);
hold on, grid on, box on
plot(parameters.time_vec, parameters.a_vec, 'LineWidth',1.5,'Color','b')
xlabel('Time [s]'), ylabel('Acceleration in the driving cycle [m/s^2]')
sp(3) = subplot(3,1,3);
hold on, grid on, box on
plot(tvec, xstar_decompressed,'LineWidth',1.5, 'LineWidth',1.5,'Color','b')
xlabel('Time [s]'), ylabel('Control variable u'), axis([-inf inf -1.2 1.2])
linkaxes(sp,'x'), clear sp

figure(5), set(gcf, 'Color', 'White'),
grid on, hold on,xlabel('time [s]'),ylabel('fuel consumption [kg]')
plot(tvec, xstar_decompressed,'LineWidth',1.5,'Color','b')
xlabel('Time [s]'), ylabel('Optimization variable x [-]'), axis([-inf inf -1.2 1.2])

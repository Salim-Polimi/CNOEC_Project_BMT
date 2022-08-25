clear all
close all
clc

%%
addpath('results')
load('50iter.mat');

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
        u_vec(ind,1) = -0.3;
    end
    
    
    [m_f_dot,soc_dot] = fuel_consumption(parameters, u_vec(ind,1),ind);
    
    m_f_dot_vec(ind,1) = m_f_dot; 
    m_f2(ind,1) = m_f2(ind-1,1)+Ts*m_f_dot; 

    soc_dot_vec(ind,1) = soc_dot;
    soc_vec2(ind,1) = soc_vec2(ind-1,1)+Ts*soc_dot;
    
    if soc_vec2(ind,1) <= 0.1
        u_vec(ind+1,1) = 0;
    end
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
figure(2), set(gcf, 'Color', 'White'),
grid on, hold on,xlabel('time [s]'),ylabel('fuel consumption [kg]')
plot(tvec, m_f1,'LineWidth',1.5,'color', 'b', 'DisplayName', 'ICE Only')
plot(tvec, m_f2,'LineWidth',1.5,'color', 'r', 'DisplayName', 'Heuristic')
plot(tvec, m_f3,'LineWidth',1.5,'color', 'g', 'DisplayName', 'Optimum Ctrl')
legend show

figure(3), set(gcf, 'Color', 'White'),
grid on, hold on, xlabel('time [s]'),ylabel('soc'), axis([-inf inf 0 1])
plot(tvec, soc_vec1,'LineWidth',1.5,'color', 'b', 'DisplayName', 'ICE Only')
plot(tvec, soc_vec2,'LineWidth',1.5,'color', 'r', 'DisplayName', 'Heuristic')
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



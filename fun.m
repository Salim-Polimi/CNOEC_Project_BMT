function [v] = fun(u_vec_compressed,parameters)

N = parameters.N;
interval_size = parameters.interval_size;
m_f(1,1) = 0;
Ts = 1;
soc_vec(1,1) = parameters.initial_soc;

% decompression of the decision variable
u_vec=[];
for i=1:length(u_vec_compressed)
    for j=1:round(N/interval_size)
        u_vec = [u_vec; u_vec_compressed(i)];
    end
end
u_vec = [u_vec; u_vec(end)];

% simulation
for ind=2:N
    [m_f_dot,soc_dot] = fuel_consumption(parameters, u_vec(ind,1),ind);    
    
    m_f_dot_vec(ind,1) = m_f_dot; 
    m_f(ind,1) = m_f(ind-1,1) + Ts*m_f_dot;
    
    soc_dot_vec(ind,1) = soc_dot;
    soc_vec(ind,1) = soc_vec(ind-1,1)+Ts*soc_dot;
end

F = [m_f(N,1)];
fun = F'*F;

% nonlinear equality constraint
g = soc_vec(N,1) - parameters.initial_soc;

% nonlinear inequality constraint h1: soc_vec(ind,1) > 0.1 for all ind
h1 = [soc_vec(1,1)-0.1];
for ind=2:N
    h1 = [h1; soc_vec(ind,1)-0.1];
end

% nonlinear inequality constraint h2: soc_vec(ind,1) < 0.9 for all ind
h2 = [0.9-soc_vec(1,1)];
for ind=2:N
    h2 = [h2; 0.9-soc_vec(ind,1)];
end

v=[fun;g;h1;h2];

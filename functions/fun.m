function [v] = fun(u_vec_compressed,parameters)

N = parameters.N;
Ts = 1;
interval_size = parameters.interval_size;
dim_decision_var = parameters.dim_decision_var;

m_f = zeros(N,1);
m_f(1,1) = 0;

soc_vec = zeros(N,1);
soc_vec(1,1) = parameters.initial_soc;

% decompression of the decision variable
u_vec = zeros(N,1);

for i=0:length(u_vec_compressed)-1
    for j=1:round(N/interval_size)
        ind = i*dim_decision_var + j;
        u_vec(ind,1) = u_vec_compressed(i+1,1);
    end
end

% simulation
for ind=2:N
    [m_f_dot,soc_dot] = fuel_consumption(parameters, u_vec(ind,1),ind);    
    
    % m_f_dot_vec(ind,1) = m_f_dot; 
    m_f(ind,1) = m_f(ind-1,1) + Ts*m_f_dot;
    
    % soc_dot_vec(ind,1) = soc_dot;
    soc_vec(ind,1) = soc_vec(ind-1,1)+Ts*soc_dot;
end

F = m_f(N,1);

% nonlinear equality constraint
g = soc_vec(N,1) - parameters.initial_soc;

% nonlinear inequality constraint h1: soc_vec(ind,1) > 0.1 for all ind
h1 = zeros(N,1);
h1(1,1) = soc_vec(1,1)-0.1;

for ind=2:N
    h1(ind,1) = soc_vec(ind,1)-0.1;
end

% nonlinear inequality constraint h2: soc_vec(ind,1) < 0.9 for all ind
h2 = zeros(N,1);
h2(1,1) = 0.9-soc_vec(1,1);

for ind=2:N
    h2(ind,1) = 0.9-soc_vec(ind,1);
end

v=[F;g;h1;h2];

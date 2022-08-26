function [m_f_dot,soc_dot] = fuel_consumption(parameters,u,ind)

% import parameters
rho_a   =   parameters.rhoa;
Af      =   parameters.Af;
Cd      =   parameters.Cd;
g       =   parameters.g;
Crr     =   parameters.Crr;
m_v     =   parameters.m_v;
R       =   parameters.R;
lambda  =   parameters.lambda;
eta_em  =   parameters.eta_em;
Voc     =   parameters.Voc;
R0      =   parameters.R0;
Qnom    =   parameters.Qnom;
eta_coul=   parameters.eta_coul;
lasso_param = parameters.lasso_param;

reg_break_limit = parameters.reg_break_limit;

vref=parameters.v_vec(ind);
a=parameters.a_vec(ind);

% Pw Computation
Finertia = m_v*a;
Froll = Crr*m_v*g;
Fdrag = 1/2*rho_a*Af*Cd*vref^2;
Fw = Finertia + Froll + Fdrag;
Pw = Fw*vref;

% Power Split Computation
if Pw >= 0
    % Traction / Costing
    Pem = Pw*u;
    Peng = Pw*(1-u);
else
    % Regenerative Braking
    if abs(Pw) < abs(reg_break_limit)
        Pem = Pw;
    else
        Pem = -abs(reg_break_limit);
    end
    Peng = 0;
end

% Fuel Consumption (m_f_dot)
if vref==0
    Pf = 0; % "start&stop not controled by optimization" hypothesis
else
    ww=vref/R;
    Tw=Peng/ww;          
    tau = gearbox(vref);
    weng=ww*tau;
    Teng=Tw/tau;
    eta_eng=0.01*efficiency_computation(weng,Teng,lasso_param);
    Pf=(Peng)/eta_eng; 
end
m_f_dot=Pf/lambda;

% Battery Consumption Calculation (soc_dot)
Pb=Pem/eta_em;
Ib=Voc/(2*R0)-sqrt((Voc/(2*R0))^2-Pb/R0);
if Ib<0
    soc_dot = -Ib*(eta_coul/Qnom);
else
    soc_dot = -Ib*(1/(Qnom*eta_coul));
end



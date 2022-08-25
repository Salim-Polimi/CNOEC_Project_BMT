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
    if Pw < -10000
        Pw = -10000;
    end
    Pem = Pw;   % We need to check if this doesn't past the limits on the power battery
    Peng = 0;
end

% Fuel Consumption (m_f_dot)
if vref==0
    Pf = 0; % "start&stop not controled by optimization" hypothesis
else
    ww=vref/R;
    Tw=Peng/ww;          % ww is already in rad/s
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
soc_dot = -Ib/Qnom*eta_coul^(sign(-Ib));


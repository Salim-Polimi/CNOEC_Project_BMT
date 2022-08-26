function [eta]=efficiency_computation(vel,torque,lasso_param)
%function which compute our efficieny based on the sum of gaussians
%computed with the LASSO.

xav=lasso_param.xav;
yav=lasso_param.yav;
omega=lasso_param.omega;
np_model=lasso_param.np_model;
theta_CVX=lasso_param.theta_CVX;

PHI=zeros(1,np_model);  % Initialize regressor matrix (Nmeasurements, each time wih np_model parameters
for ind=1:np_model
    PHI(1,ind)=exp(-((vel-xav(ind)).^2 + (torque-yav(ind)).^2)/(2*omega^2));  %build the polynomial structure, starting from the largest exponential 
end
eta_map= theta_CVX'*PHI';
if eta_map>=2 %ho aggiunto questo if per evitare le singolarità sull'efficienza 
    eta=eta_map;
else
    eta=2;
end

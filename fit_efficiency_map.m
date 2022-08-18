function [xav,yav,omega,np_model,theta_CVX] = fit_efficiency_map(speed_data,torque_data,eta_data)

%discretizzazione uniforme dello spazio 
x= linspace(-0,550,1000);
y=linspace(-0,250,1000);
[X,Y]= meshgrid(x,y);   

eta_meas= eta_data'; %Nx1 column vector 
N=length(speed_data); %number of sampled data, NB: li devo ordinare in ordine crescente? 

%COUPLED (V,T) VECTOR
zvals=zeros(N,2);
for i=1:N
    zvals(i,1)= speed_data(1,i);
    zvals(i,2)= torque_data(1,i);
end

%PARAMETERS OF THE GAUSSIANS (by now I used the same sampled data)
for i=1:N
    xav(i)= speed_data(1,i);
    yav(i)= torque_data(1,i);
end

epsilon= 3; %desired bound of the lasso
np_model=220;
omega=35;

PHI=zeros(N,np_model);  % Initialize regressor matrix (Nmeasurements, each time wih np_model parameters
for ind=1:np_model
    PHI(:,ind)=exp(-((zvals(:,1)-xav(ind)).^2 + (zvals(:,2)-yav(ind)).^2)/(2*omega^2));  %build the polynomial structure, starting from the largest exponential 
end

cvx_begin  
variable theta_CVX(np_model,1) %la variabile da trovare è il valore dei parametri 
minimize norm(theta_CVX,1) 
subject to %constraints
norm(eta_meas-PHI*theta_CVX,inf)<= epsilon %maximum displacement from the real function
cvx_end %end of the min problem which leads to 'theta_CVX', much more similar to the true value 

eta_CVX               =   PHI*theta_CVX; %molto simili a quelli campionati
figure(2),plot3(zvals(:,1),zvals(:,2), eta_meas, '*'), grid on, hold on;

j=22;  %number of sampled torque for each velocities
figure(2),plot3(zvals(1:j,1),zvals(1:j,2),eta_CVX((1:j),1),'y','linewidth',2),hold on;
for i=2:(N/22) 
    figure(2),plot3(zvals(((i-1)*j+1):(j*i),1),zvals(((i-1)*j+1):(j*i),2),eta_CVX(((i-1)*j+1):(j*i)),'y','linewidth',2),hold on;
end
title('LASSO estimated points');

%PLOT OF THE FUNCTION OBTAINED WITH LASSO
fun=theta_CVX(1)*exp(-((X-xav(1)).^2 + (Y-yav(1)).^2)/(2*omega^2));
for ind=2:np_model
    fun=fun+theta_CVX(ind)*exp(-((X-xav(ind)).^2 + (Y-yav(ind)).^2)/(2*omega^2));
end
  
Z_CVX=fun; 
figure(3)
surf(X,Y,Z_CVX)
title('LASSO curve fitting');
shading interp
axis tight

end
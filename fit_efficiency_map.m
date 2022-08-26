clear all 
close all
clc

addpath('fit_efficiency_map')

% GOAL:  create a  continuous function which associates each operating 
%        condition of the ICE to an efficiency.
%        Fit the sampled points with a sum of 3D gaussian functions.

% 1.a)   LS FITTING 
% 1.b)   LS FITTING with 110 parameters to simulate overfitting results
% 2.a)   LASSO optimization with all 220 points. 
% 2.b)   LASSO optimization with 110 points. Other 110 left for validation.


% Import data from a given efficiency map
load('sampled_data.mat')
                                                  
speed_data    =      sampled_data(1,:);                                    %engine speed (RPM)
torque_data   =      sampled_data(2,:);                                    %engine torque(Nm)
eta_data      =      sampled_data(3,:);                                    %efficiencies, eta=f(speed, torque)

eta_meas      =      eta_data';                                            %Nx1 column vector 
M             =      length(speed_data);                                   %number of sampled data


%BUILD (w_eng, T_eng) VECTOR
zvals          =      zeros(M,2);
for i=1:M
    zvals(i,1) =      speed_data(1,i);
    zvals(i,2) =      torque_data(1,i);
end


%% 1.a) LS with 110 points and 55 parameters
%  Use 110 points for fitting and other 110 for the validation
%  Use only 55 parameters in order to avoid overfitting


np_model      =     55;                                                    %number of model parameters
M             =     110;                                                   %number of used points 

%averages of the gaussian functions
xav1    =  zeros(1,110);
yav1    =  zeros(1,110);
for i=1:length(speed_data)
    if mod(i,4)==0
       xav1(i/4) =      speed_data(1,i);
       yav1(i/4) =      torque_data(1,i);
    end
end

%divide into FITTING SET and VALIDATION SET         
%BUILD VECTOR (w_eng, T_eng) AND EFFICIENCIES FOR LS

zvals_validation     = [];
eta_meas_validation  = [];
zvals_fitting        = zeros(110,2);
eta_meas_fitting     = zeros(110,1);

for i=1:length(zvals)
    if mod(i,2)==0
       zvals_fitting(i,:)       =   zvals(i,:);
       eta_meas_fitting(i,1)    =   eta_meas(i);
    else
       zvals_validation         =   [zvals_validation;zvals(i,:)];
       eta_meas_validation      =   [eta_meas_validation;eta_meas(i)];
    end
end

zvals_fitting      =    zvals_fitting(~all(zvals_fitting == 0, 2),:);
eta_meas_fitting   =    eta_meas_fitting(~all(eta_meas_fitting == 0, 2),:);

omega  =   30;          %variance of the gaussian functions 

%BUILD THE REGRESSOR MATRIX
PHI=zeros(M,np_model);                                                     
for ind=1:np_model
    PHI(:,ind)  = exp(-((zvals_fitting(:,1)-xav1(ind)).^2 + (zvals_fitting(:,2)-yav1(ind)).^2)/(2*omega^2));  %build the polynomial structure, starting from the largest exponential 
end

A_LS            =   pinv(PHI)*eta_meas_fitting;    

%plot all points etimated efficiencies
PHI=zeros(220,np_model);                                                     
for ind=1:np_model
    PHI(:,ind)  = exp(-((zvals(:,1)-xav1(ind)).^2 + (zvals(:,2)-yav1(ind)).^2)/(2*omega^2));  %build the polynomial structure, starting from the largest exponential 
end
eta_LS     =   PHI*A_LS; 

figure(10),plot3(zvals(:,1),zvals(:,2), eta_meas, '*'), grid on, hold on;
j          =   22;                   
figure(10),plot3(zvals(1:j,1),zvals(1:j,2),eta_LS((1:j),1),'*', 'Color', 'g'),hold on;

for i=2:(220/22) 
figure(10),plot3(zvals(((i-1)*j+1):(j*i),1),zvals(((i-1)*j+1):(j*i),2),eta_LS(((i-1)*j+1):(j*i)),'*', 'Color', 'g'),hold on;
end
legend('sampled points', 'estimated function')
xlabel('engine angular speed(RPM*10)')
ylabel('engine torque(Nm)')
zlabel('engine efficiency (%)')
title('Least Squares estimated points');


%SPACE DISCRETIZATION
x               =   linspace(-0,550,1000);
y               =   linspace(-0,250,1000);
[X,Y]           =   meshgrid(x,y);
 

%EFFICIENCY COMPUTATION FOR DISCRETIZED (w_eng, T_eng) POINTS 
fun=A_LS(1)*exp(-((X-xav1(1)).^2 + (Y-yav1(1)).^2)/(2*omega^2));
  for ind=2:np_model
       fun    =   fun+A_LS(ind)*exp(-((X-xav1(ind)).^2 + (Y-yav1(ind)).^2)/(2*omega^2));
  end
  
Z_LS          =   fun; 

figure
surf(X,Y,Z_LS)
title('Least Squares efficiency map');
xlabel('engine angular speed(RPM*10)')
ylabel('engine torque(Nm)')
zlabel('engine efficiency(%)')
shading interp
axis tight
 
%VALIDATION
X                        =   zvals_validation(:,1);
Y                        =   zvals_validation(:,2);


fun_validation=A_LS(1)*exp(-((X-xav1(1)).^2 + (Y-yav1(1)).^2)/(2*omega^2));
  for ind=2:np_model
      fun_validation     =   fun_validation+A_LS(ind)*exp(-((X-xav1(ind)).^2 + (Y-yav1(ind)).^2)/(2*omega^2));
  end

err_LS                   =   sqrt((eta_meas_validation-fun_validation)'*(eta_meas_validation-fun_validation))


%% 1.b) Simulated overfitting (LS with 110 points and 110 parameters)
%  example of how the LS estimate suffers of overfitting if the dimension 
%  of the parameters vector is doubled.

 np_model      =      110;                                                 %number of model parameters   
 M             =      110;                                                 %number of used points
 
for i=1:length(speed_data)
    if mod(i,2)==0
    xav1(i/2)  =      speed_data(1,i);
    yav1(i/2)  =      torque_data(1,i);
    end
end

zvals_fitting      =      zeros(110,2);
eta_meas_fitting   =      zeros(110,1);
for i=1:length(zvals)
    if mod(i,2)==0
       zvals_fitting(i,:)    =      zvals(i,:);
       eta_meas_fitting(i,1) =      eta_meas(i);
    end
end
zvals_fitting       =    zvals_fitting(~all(zvals_fitting == 0, 2),:);
eta_meas_fitting    =    eta_meas_fitting(~all(eta_meas_fitting == 0, 2),:);

omega         =       30;                                                  %variance of gaussians function

PHI=zeros(M,np_model);                                                     %Initialize regressor matrix 
for ind=1:np_model
    PHI(:,ind) = exp(-((zvals_fitting(:,1)-xav1(ind)).^2 + (zvals_fitting(:,2)-yav1(ind)).^2)/(2*omega^2));  %build the polynomial structure, starting from the largest exponential 
end



A_LS          =   pinv(PHI)*eta_meas_fitting;    
eta_LS        =   PHI*A_LS;

%SPACE DISCRETIZATION
x             =   linspace(-0,550,1000);
y             =   linspace(-0,250,1000);
[X,Y]         =   meshgrid(x,y);
 


%EFFICIENCY COMPUTATION FOR DISCRETIZED (w_eng, T_eng) POINTS 
fun=A_LS(1)*exp(-((X-xav1(1)).^2 + (Y-yav1(1)).^2)/(2*omega^2));
  for ind=2:np_model
       fun    =   fun+A_LS(ind)*exp(-((X-xav1(ind)).^2 + (Y-yav1(ind)).^2)/(2*omega^2));
   end
  
Z_LS          =   fun; 

figure
surf(X,Y,Z_LS)
title('Least Squares efficiency map');
xlabel('engine angular speed(RPM*10)')
ylabel('engine torque(Nm)')
zlabel('engine efficiency(%)')
shading interp
axis tight


%% 2.a) LASSO with 220 points
epsilon      =    3;                                                       %maximum fitting error
np_model     =    220;                                                     %dimension of the parameter vector
omega        =    35;                                                      %gaussian variance 
M            =    220;                                                     %number of used points
PHI          =    zeros(M,np_model);                                       %Initialize regressor matrix (Nmeasurements, each time wih np_model parameters

xav          =    zeros(1,220);
yav          =    zeros(1,220);
for i=1:length(speed_data)
    xav(i)    =      speed_data(1,i);
    yav(i)    =      torque_data(1,i);
end

%BUILD THE POLINOMIAL STRUCTURE 
for ind=1:np_model
    PHI(:,ind)=exp(-((zvals(:,1)-xav(ind)).^2 + (zvals(:,2)-yav(ind)).^2)/(2*omega^2));  
end

cvx_begin  
variable theta_CVX(np_model,1)               %theta_CVX= A_LASSO in the report    
minimize norm(theta_CVX,1) 
%constraints
subject to                            
norm(eta_meas-PHI*theta_CVX,inf)<= epsilon;  
cvx_end 

eta_CVX      =   PHI*theta_CVX; 

%PLOT OF THE FUNCTION OBTAINED WITH LASSO
x             =   linspace(-0,550,1000);
y             =   linspace(-0,250,1000);
[X,Y]         =   meshgrid(x,y);

fun = theta_CVX(1)*exp(-((X-xav(1)).^2 + (Y-yav(1)).^2)/(2*omega^2));
  for ind=2:np_model
       fun=fun+theta_CVX(ind)*exp(-((X-xav(ind)).^2 + (Y-yav(ind)).^2)/(2*omega^2));
   end
  
Z_CVX   =  fun; 

figure
surf(X,Y,Z_CVX)
title('LASSO eta(weng,Teng)');
shading interp
axis tight



%% LASSO with 110 points (other 110 for validation)
epsilon      =    3;                  %desired bound of the lasso
np_model     =    110;
omega        =    35;
M=110;

%divide points into fitting and validation set
for i=1:length(speed_data)
    if mod(i,2)==0
       xav1(i/2) =      speed_data(1,i);
       yav1(i/2) =      torque_data(1,i);
    end
end

zvals_validation         = [];
eta_meas_validation      = [];
zvals_fitting            = zeros(110,2);
eta_meas_fitting         = zeros(110,1);

for i=1:length(zvals)
    if mod(i,2)==0
       zvals_fitting(i,:)      =      zvals(i,:);
       eta_meas_fitting(i,1)   =      eta_meas(i);
    else
       zvals_validation        =      [zvals_validation;zvals(i,:)];
       eta_meas_validation     =      [eta_meas_validation;eta_meas(i)];
    end
end

zvals_fitting       =    zvals_fitting(~all(zvals_fitting == 0, 2),:);
eta_meas_fitting    =    eta_meas_fitting(~all(eta_meas_fitting == 0, 2),:);


M               =    110;

PHI             =    zeros(M,np_model); 

%BUILD THE POLINOMIAL STRUCTURE 
for ind=1:np_model
    PHI(:,ind)=exp(-((zvals_fitting(:,1)-xav1(ind)).^2 + (zvals_fitting(:,2)-yav1(ind)).^2)/(2*omega^2));  
end

cvx_begin  
variable theta_CVX(np_model,1)        
minimize norm(theta_CVX,1) 

subject to                            
norm(eta_meas_fitting-PHI*theta_CVX,inf)<= epsilon;  
cvx_end 

eta_CVX      =   PHI*theta_CVX; 

%plot all points estimated efficiency

PHI             =    zeros(220,np_model); 
for ind=1:np_model
    PHI(:,ind)=exp(-((zvals(:,1)-xav1(ind)).^2 + (zvals(:,2)-yav1(ind)).^2)/(2*omega^2));  
end
eta_CVX      =   PHI*theta_CVX; 

figure(3),plot3(zvals(:,1),zvals(:,2), eta_meas, '*'), grid on, hold on;
j            =   22;                   
figure(3),plot3(zvals(1:j,1),zvals(1:j,2),eta_CVX((1:j),1),'*', 'Color', 'r'),hold on;

for i=2:(220/22) 
figure(3),plot3(zvals(((i-1)*j+1):(j*i),1),zvals(((i-1)*j+1):(j*i),2),eta_CVX(((i-1)*j+1):(j*i)),'*', 'Color', 'r'),hold on;
end
legend('sampled points', 'estimated function')
xlabel('engine angular speed(RPM*10)')
ylabel('engine torque(Nm)')
zlabel('engine efficiency(%)')
title('LASSO estimated points');


%PLOT OF THE FUNCTION OBTAINED WITH LASSO
x             =   linspace(-0,550,1000);
y             =   linspace(-0,250,1000);
[X,Y]         =   meshgrid(x,y);

fun = theta_CVX(1)*exp(-((X-xav1(1)).^2 + (Y-yav1(1)).^2)/(2*omega^2));
  for ind=2:np_model
       fun=fun+theta_CVX(ind)*exp(-((X-xav1(ind)).^2 + (Y-yav1(ind)).^2)/(2*omega^2));
   end
  
Z_CVX   =  fun; 

figure
surf(X,Y,Z_CVX)
title('LASSO efficiency map');
xlabel('engine angular speed(RPM*10)')
ylabel('engine torque(Nm)')
zlabel('engine efficiency(%)')
shading interp
axis tight


%VALIDATION

X     =   zvals_validation(:,1);
Y     =   zvals_validation(:,2);


fun_validation = theta_CVX(1)*exp(-((X-xav1(1)).^2 + (Y-yav1(1)).^2)/(2*omega^2));
  for ind=2:np_model
      fun_validation=fun_validation+theta_CVX(ind)*exp(-((X-xav1(ind)).^2 + (Y-yav1(ind)).^2)/(2*omega^2));
  end

err_LASSO       =   sqrt((eta_meas_validation-fun_validation)'*(eta_meas_validation-fun_validation))







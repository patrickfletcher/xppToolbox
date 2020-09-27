% script demonstrating usage of ode2m and ode_euler to solve a system with
% Wiener variables

clear

odefilename='lactotroph_noise.ode';

[mFunctionName, xppdata]=ode2m(odefilename);

%lists of parameter, variable, and aux variable names:
parNames=xppdata.parNames; 
varNames=xppdata.varNames;
auxNames=xppdata.auxNames;

p=xppdata.p0; %the default parameter vector

% one way to change parameter values is like so:
p(parNames=="gbk")=0.0;

% initial condition
y0=xppdata.x0;

% access XPP options that were found in the ODE file
dt=xppdata.opt.dt; 
total=xppdata.opt.total;

%% call matlab solver
%define the RHS anonymous function in matlab - must be done *after* changing any parameters of interest
odefun=eval(['@(t,y,w) ' mFunctionName '(t,y,w,p)']);
[t,y,dy,aux,w]=ode_euler(odefun,[0,total],dt,y0,xppdata.nWiener);
%w is the sequence of random normals used in the simulation

% plot result
figure(1)
subplot(2,1,1)
var=varNames=="v";
plot(t,y(:,var))
xlabel('t')
ylabel(varNames(var))

subplot(2,1,2)
auxvar=auxNames=="ica";
plot(t,aux(:,auxvar))
xlabel('t')
ylabel(auxNames(auxvar))
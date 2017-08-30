% script to demonstrate xppConverter ode to m-file usage
clear

odefilename='lactotroph.ode';

[mFunctionName, xppdata]=xppConverter('m',odefilename);

%parameter values can be changed like so:
p=xppdata.p0; %first get the default parameter vector
p(strcmp(xppdata.parNames,'gbk'))=0.25; %change a value of interest

%define the RHS anonymous function in matlab
odefun=eval(['@(t,y) ' mFunctionName '(t,y,p)']);

%initial condition
y0=xppdata.x0;

total=xppdata.opt.total;

%call matlab solver
[t,y]=ode45(odefun,[0,total],y0);

% plot result
varIx=1;

figure
plot(t,y(:,varIx))
xlabel('t')
ylabel(xppdata.varNames(varIx))

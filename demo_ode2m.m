% script demonstrating basic usage of ode2m on an ODE file
clear

odefilename='lactotroph.ode';

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
%define the RHS anonymous function in matlab - must be done *after*
%changing any parameters of interest. The anonymous function will store the
%values of p as they are once reaching this line of code.
odefun=eval(['@(t,y) ' mFunctionName '(t,y,p)']);
[t,y]=ode45(odefun,[0,total],y0);

%compute the auxiliary variables from the solution values
aux=zeros(length(t),xppdata.nAux);
for i=1:length(t)
    [~,tmp]=odefun(t(i),y(i,:));
    aux(i,:)=tmp';
end

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
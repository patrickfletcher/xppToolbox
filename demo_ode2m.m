% script to demonstrate xppConverter ode to m-file usage
clear

odefilename='lactotroph.ode';

[mFunctionName, xppdata]=ode2m(odefilename);

%parameter values can be changed like so:
p=xppdata.p0; %first get the default parameter vector
parNames=xppdata.parNames;

p(parNames=="gbk")=0.0; %change a value of interest

%initial condition
y0=xppdata.x0;
varNames=xppdata.varNames;

%access XPP options that were found in the ODE file
dt=xppdata.opt.dt; 
total=xppdata.opt.total;

%% call matlab solver
if xppdata.nWiener==0
    
    %define the RHS anonymous function in matlab - must be done *after* changing any parameters of interest
    odefun=eval(['@(t,y) ' mFunctionName '(t,y,p)']);
    [t,y]=ode45(odefun,[0,total],y0);
    
    %compute the auxiliary variables (if not using ode_euler)
    aux=zeros(length(t),xppdata.nAux);
    for i=1:length(t)
        [~,tmp]=odefun(t(i),y(i,:));
        aux(i,:)=tmp';
    end

else
    odefun=eval(['@(t,y,w) ' mFunctionName '(t,y,w,p)']);
    [t,y,dy,aux,w]=ode_euler(odefun,[0,total],dt,y0,xppdata.nWiener);
    %w is the sequence of random normals used in the simulation
end


% plot result
figure(1)
var=varNames=="v";
plot(t,y(:,var))
xlabel('t')
ylabel(varNames(var))

figure(2)
auxNames=xppdata.auxNames;
plot(t,aux(:,1))
xlabel('t')
ylabel(auxNames(1))
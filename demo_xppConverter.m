% script to demonstrate xppConverter ode to m-file usage
clear

odefilename='lactotroph.ode';

[mFunctionName, xppdata]=xppConverter(odefilename);

%parameter values can be changed like so:
p=xppdata.p0; %first get the default parameter vector
p(strcmp(xppdata.parNames,'gbk'))=0.25; %change a value of interest


%initial condition
y0=xppdata.x0;

total=xppdata.opt.total;

%% call matlab solver
if xppdata.nWiener==0
    %define the RHS anonymous function in matlab
    odefun=eval(['@(t,y) ' mFunctionName '(t,y,p)']);
    [t,y]=ode45(odefun,[0,total],y0);
else
    odefun=eval(['@(t,y,w) ' mFunctionName '(t,y,w,p)']);
    [t,y]=ode_euler(odefun,[0,total],0.5,y0,xppdata.nWiener);
end

% plot result
varIx=1;

figure(1)
plot(t,y(:,varIx))
xlabel('t')
ylabel(xppdata.varNames(varIx))

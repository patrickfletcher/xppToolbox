function [t,Y,DY,AUX,w]=ode_euler(getRHS,tspan,dt,y0,nWiener)
%ODE_EULER A simple euler solver supporting Wiener variables.
%
% 

ti=tspan(1);
yi=y0(:);
nSteps=ceil( (tspan(2)-tspan(1))/dt )+1;
t=zeros(nSteps,1);
Y=zeros(nSteps,length(y0));
DY=zeros(nSteps,length(y0));
w=zeros(nSteps,length(y0));

%initial values
t(1)=ti;
Y(1,:)=yi(:)';

w(1,:)=randn(nWiener,1)/sqrt(dt);

[dyi, tmpaux]=getRHS(ti, Y(1,:), w(1,:));
DY(1,:)=dyi(:)';

AUX=zeros(nSteps,length(tmpaux));
AUX(1,:)=tmpaux(:)';

%timestep loop
for i=2:nSteps
    
    w(i,:)=randn(nWiener,1)/sqrt(dt);
    [k1, aux]=getRHS(ti, yi, w(i,:)); 
    yi=yi+k1*dt;
    ti=t(1)+(i-1)*dt;
    
    %store
    t(i,1)=ti;
    Y(i,:)=yi(:)'; %(:)' forces row vector
    DY(i,:)=dyi(:)';
    AUX(i,:)=aux(:)';
end
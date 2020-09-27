function [t,Y,DY,AUX,W]=ode_euler(getRHS,tspan,dt,y0,nWiener)
%ODE_EULER A simple euler solver supporting Wiener variables.

%initial values
t=tspan(1):dt:tspan(2);
nSteps=length(t)-1;
Y=zeros(nSteps+1,length(y0));
Y(1,:)=y0(:)';

w=randn(nWiener,1)/sqrt(dt);
[dyi, tmpaux]=getRHS(t(1), Y(1,:), w);

if nargout>2
DY=zeros(nSteps+1,length(y0));
W=zeros(nSteps+1,length(y0));
AUX=zeros(nSteps+1,length(tmpaux));

DY(1,:)=dyi(:)';
AUX(1,:)=tmpaux(:)';
W(1,:)=w;
end

%timestep loop
yi=y0(:);
for i=2:nSteps+1
    
    w=randn(nWiener,1)/sqrt(dt);
    [k1, aux]=getRHS(t(i), yi, w); 
    yi=yi+k1*dt;
    Y(i,:)=yi(:)'; %(:)' forces row vector
    
    if nargout>2
        DY(i,:)=dyi(:)';
        AUX(i,:)=aux(:)';
        W(i,:)=w;
    end
    
    if mod(i,round(nSteps/10))==0
        fprintf('.')
    end
end
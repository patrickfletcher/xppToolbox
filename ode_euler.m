function [t,Y,DY,AUX,W]=ode_euler(getRHS,tspan,dt,y0,nWiener)
%ODE_EULER A simple euler solver supporting Wiener variables.
%
% 

% ti=tspan(1);
yi=y0(:);
nSteps=ceil( (tspan(2)-tspan(1))/dt )+1;
% t=zeros(nSteps,1);
Y=zeros(nSteps,length(y0));

%initial values
% t(1)=ti;
t=tspan(1):dt:nSteps*dt;
Y(1,:)=yi(:)';


w=randn(nWiener,1)/sqrt(dt);
[dyi, tmpaux]=getRHS(t(1), Y(1,:), w);

if nargout>2
DY=zeros(nSteps,length(y0));
W=zeros(nSteps,length(y0));
AUX=zeros(nSteps,length(tmpaux));

DY(1,:)=dyi(:)';
AUX(1,:)=tmpaux(:)';
W(1,:)=w;
end

%timestep loop
for i=2:nSteps+1
    
    w=randn(nWiener,1)/sqrt(dt);
    [k1, aux]=getRHS(t(i), yi, w); 
    yi=yi+k1*dt;
%     ti=t(1)+(i-1)*dt;
    
    %store
%     t(i,1)=ti;
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
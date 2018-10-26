function result=heaviside(x,useHalfMax)
result=0;
if ~exist('useHalfMax','var'), useHalfMax=false; end
if x>0
    result=1;
elseif x==0 && useHalfMax
    result=0.5;
end
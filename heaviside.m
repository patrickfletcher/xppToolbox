function result=heaviside(x, useHalfMax)
arguments
    x
    useHalfMax=false
end
result=zeros(size(x));
result(x>0)=1;
if useHalfMax
    result(x==0)=0.5;
end
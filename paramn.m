function [output,a,b]=paramn(input)
%%%parameterization with normalization
% [m,n]=size(input);
% temp=reshape(input',m*n,1);
temp=input/norm(input);
a=temp(1);
b=temp(2:end);
output=2/csinc(acos(a))*b;
if norm(output)>pi
    output=(1-2*pi/norm(output)*ceil((norm(output)-pi)/(2*pi)))*output;
end
end
function out=derivaparam(in)
[inhat,~,b]=paramn(in);
out21=-0.5*b';
if (norm(inhat)==0)
    out22=1/2*eye(size(b,1));
else
    out22=csinc(norm(inhat)/2)/2*eye(size(b,1))+1/(4*norm(inhat))*dcsinc(norm(inhat)/2)...
        *inhat*inhat';
end
out=[out21;out22];


function out=dehomo(in)
n=size(in,1);
switch n
    case 3;
        out=[in(1,:)./in(3,:);
            in(2,:)./in(3,:)];
    case 4;
        out=[in(1,:)./in(4,:);
            in(2,:)./in(4,:);
            in(3,:)./in(4,:)];
end

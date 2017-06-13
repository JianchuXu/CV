function [ imJ ] = nonmax( iminput, wsize )
temp = (wsize-1) /2;
temp111 = padarray(iminput,[temp,temp]);
padding = zeros(size(temp111));

for row = (1+temp):(size(iminput,1)+temp)
    for col = (1+temp):(size(iminput,2)+temp)
        if temp111(row,col) == max(max(temp111(row-temp:row+temp,col-temp:col+temp)))
           padding(row,col) = 1;
        end
    end
end

imJ = padding((1+temp):(size(iminput,1)+temp),...
    (1+temp):(size(iminput,2)+temp)) .*...
    temp111((1+temp):(size(iminput,1)+temp),...
    (1+temp):(size(iminput,2)+temp)) ;
end


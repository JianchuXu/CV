function [ minorI ] = minor( im1, im2, wsize ,thres)
% compute the minor image.
im1 = double(im1);im2 = double(im2);
h = ones(wsize,wsize);
minor1 = conv2(double(im1 .* im1),double(h),'same')/(wsize^2);
minor2 = conv2(double(im2 .* im2),double(h),'same')/(wsize^2);
minor3 = conv2(double(im1 .* im2),double(h),'same')/(wsize^2);
minorI = zeros(size(im1));
for i = 1:size(im1,1)
    for j = 1:size(im1,2)
        N = double([minor1(i,j) minor3(i,j); minor3(i,j) minor2(i,j)]);
        lambdamin = (trace(N) - sqrt(power(trace(N),2) - 4*det(N)))/2;
        if lambdamin > thres
            minorI(i,j) = lambdamin;
        end
    end
end
end


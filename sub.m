function [ xs, ys] = sub( im1, im2, x, y, wsize )
im1 = double(im1);im2 = double(im2);
b1 = zeros(size(im1));b2 = zeros(size(im1));
for i = 1:size(im1,1)
    for j = 1:size(im1,2)
        b1(i,j) = im1(i,j) * im1(i,j) * i + im1(i,j) * im2(i,j) * j;
        b2(i,j) = im2(i,j) * im2(i,j) * j + im1(i,j) * im2(i,j) * i;
    end
end
h = ones(wsize);
minor1 = conv2(double(im1 .* im1),double(h),'same');
minor2 = conv2(double(im2 .* im2),double(h),'same');
minor3 = conv2(double(im1 .* im2),double(h),'same');
b1 = conv2(double(b1),double(h),'same');
b2 = conv2(double(b2),double(h),'same');

xs = zeros(size(x));
ys = zeros(size(y));
for i = 1:size(x,1)
    B = [b1(x(i),y(i));b2(x(i),y(i))];
    A = [minor1(x(i),y(i)) minor3(x(i),y(i)); 
        minor3(x(i),y(i)) minor2(x(i),y(i))];
    result = A\B;
    xs(i) = result(1);
    ys(i) = result(2);
end

end


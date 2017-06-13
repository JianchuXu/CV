 for i=1:6
subplot(2,3,i);
plot(range,errorb(i*11-10:i*11));
xlim([1 64]);
end
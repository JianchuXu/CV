function result=dcsinc(input)
if input==0
   result=0;
else
    result=cos(input)/input-sin(input)/(input^2);
end
end
function entropy = SampEn(data,r)
% 
l = length(data);
Nn = 0;
Nd = 0;
for i = 1:l-2
    for j = i+1:l-2
        if abs(data(i)-data(j))<r && abs(data(i+1)-data(j+1))<r
            Nn = Nn+1;
            if abs(data(i+2)-data(j+2))<r
                Nd = Nd+1;
            end
        end
    end
end
entropy = -log(Nd/Nn);
end
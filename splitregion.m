%Copyright (c) 2020 
%3D Systems Packaging Research Center (PRC), Georgia Tech.
%Split D-dimensional input sample space to 2^D regions.
function output_domains = splitregion(domain,d)
n_region = 2^d;

output_domains = zeros(d,2,n_region);
temp = zeros(d,3);
for i = 1:d
    temp(i,3) = domain(i,2);
    temp(i,1)  = domain(i,1);
    temp(i,2) = (temp(i,3)+temp(i,1))/2;
end

i = 1;
count = 0;
while (count < n_region)
%     index = decimalToBinaryVector(count,d);
    index = fliplr(de2bi(count,d));
    index = index + 1;
    index = flip(index);
    for a = 1:d
        output_domains(a,:,count+1) = [temp(a,index(a)) temp(a,index(a)+1)]; 
    end
    count = count + 1;
end

end

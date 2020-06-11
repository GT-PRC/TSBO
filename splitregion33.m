%Copyright (c) 2020 
%3D Systems Packaging Research Center (PRC), Georgia Tech.
%Split D-dimensional input sample space to 3 regions along the largest parameter.
function output_domains = splitregion33(domain,d)
n_region = 3;

output_domains = zeros(d,2,n_region);
temp = zeros(d,4);
spaces = domain(:,2) - domain(:,1);
for i = 1:d
    difff = (domain(i,2)-domain(i,1))/3;
    temp(i,1)  = domain(i,1);
    temp(i,2)  = domain(i,1) + difff;
    temp(i,3)  = domain(i,1) + 2*difff;
    temp(i,4)  = domain(i,1) + 3*difff;
end
[~,longest_i] = max(spaces);

for i = 1:n_region
    for a = 1:d
        if(a ~= longest_i)
            output_domains(a,:,i) = [domain(a,1) domain(a,2)];
        else
            output_domains(a,:,1) = [temp(longest_i,1) temp(longest_i,2)];
            output_domains(a,:,2) = [temp(longest_i,2) temp(longest_i,3)];
            output_domains(a,:,3) = [temp(longest_i,3) temp(longest_i,4)];

        end
    end
end
end

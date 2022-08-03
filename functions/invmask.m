function [imask] = invmask(mask)

[m,n] = size(mask);

for i=1:m
    
    for j=1:n
        if mask(i,j) == true
            mask(i,j) = false;
        else
            mask(i,j) = true;
        end
    end
    
end
imask = mask;
end
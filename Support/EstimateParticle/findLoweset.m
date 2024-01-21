function [test, index] = findLoweset(E1)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

index = zeros(size(E1,1),size(E1,2));

for i = 1 : size(E1,1)
    for j = 1 : size(E1,2)
       [test, index_t]  = min(E1(i,j,:));       
       if test == 0
           index(i,j) = 0;
       else
           index(i,j) = index_t;
       end    
    end
end


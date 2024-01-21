function weight = WeightFunc(v, x,r,index_control)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

m_v = size(v,1);
m_x = size(x,1); 

weight = zeros(m_v, m_x);
for i = 1 : m_v
 for j = 1 : m_x
    index = find(index_control==i);
    
    if(index)
   
         weight(i,index) = 1.0;
         continue;
     else

         d_2 = norm(v(i,:)-x(j,:))^2/ r(i)^2;
           %d_2 = norm(v(i,:)-x(j,:))^2/ r(j)^2;
         weight(i,j) = max(0, (1-d_2)^3);
      end 
    
 end
end


div_tmp = sum(weight,2);

 for i = 1 : m_v
     for j = 1 : m_x
            if div_tmp(i) == 0
               [idx,Distance] = knnsearch(x,v(i,:),'k',3);
                weight(i,idx(1)) = 0.5;
                weight(i,idx(2)) = 0.3;
                weight(i,idx(3)) = 0.2;
                div_tmp(i)=1;
                disp('here');
            end
            weight(i,j) = weight(i,j)/div_tmp(i);
     end
 end

end


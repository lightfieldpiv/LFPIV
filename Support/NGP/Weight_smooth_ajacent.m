function weight = Weight_smooth_ajacent(  control_point,r )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%value_test = zeros(control_point,control_point);

m_v = size(control_point,1);
m_x = size(control_point,1); 
%for i = 1 : m

%i = 1 : m_v;
value_tmp = zeros(m_x,m_v);
value = zeros(m_x,m_v);


weight = zeros(m_v,m_x);
%%
for i = 1 : m_v
 for j = 1 : m_x
     if(i~=j)
        d_2 = norm(control_point(i,:)-control_point(j,:))^2/(r(i)+r(j))^2;
        weight(i,j) = max(0,(1-d_2)^3);
     end
    
 end
end

div_tmp = sum(weight,2);

 for i = 1 : m_v
     for j = 1 : m_x
            if div_tmp(i) == 0
               [idx,Distance] = knnsearch(control_point,control_point(i,:),'k',3);
                weight(i,idx(1)) = 0.5;
                weight(i,idx(2)) = 0.3;
                weight(i,idx(3)) = 0.2;
                div_tmp(i)=1;
                disp('overall weight here');
            end
            weight(i,j) = weight(i,j)/div_tmp(i);
     end
 end



end


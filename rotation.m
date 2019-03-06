%ADV×ø±ê±ä»»
function [U,V,W]=rotation(U,V,W)
B=atan(V./U);
VS=(V.^2+U.^2).^0.5;
Q=-atan(W./VS);
Us=[U V W];
Ur=zeros(length(U),3);
for i=1:length(U)
Ur(i,1:3)=Us(i,1:3)*...
    [cos(Q(i)).*cos(B(i)) cos(Q(i)).*sin(B(i)) sin(Q(i));...
    -sin(B(i)) cos(B(i)) 0;...
    -sin(Q(i)).*cos(B(i)) -sin(Q(i)).*sin(B(i)) cos(Q(i))];
end
U=Ur(:,1);
V=Ur(:,2);
W=Ur(:,3);
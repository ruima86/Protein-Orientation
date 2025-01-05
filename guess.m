function g = guess(u,r0,p,rb)

% dy = zeros(9,1);    
% psi = y(1);   
% dpsi = y(2);
% r = y(3);
% z = y(4);
% alpha = y(5);
% beta = y(6);
% h = y(7);
% gamma = y(8);
% area = y(9);

% geometry of a sphere

g = [pi*u;...
     pi;...
     r0+rb.*sin(pi.*u);...
     rb+rb.*cos(pi.*u);...
     1;...
     1;...
     pi.*rb;...
     -p*rb/2;...
     2.*rb.*(r0.*u-rb.*(pi)^(-1).*cos(u.*pi)+rb/pi)];
 
end
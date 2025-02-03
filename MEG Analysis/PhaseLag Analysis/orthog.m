function [o,R,B]=orthog(x,y)

% RH 1999

% orthog x wrt y

if(size(x,1)==1)
	x=x';
end
if(size(y,1)==1)
	y=y';
end

B = pinv(y)*x;
R = eye(size(y,1)) - y*pinv(y);
o = R*x;
	

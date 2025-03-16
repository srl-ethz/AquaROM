function T = Tu3_conf3(rho,x01,x02,x03,x04,x05,x06)
T(:,:,1) = [];
T(:,:,2) = [];

T(:,:,3) = zeros(6,6);
T(:,:,4) = zeros(6,6);
T(:,:,5) = T(:,:,1);
T(:,:,6) = T(:,:,2);

T = 0.5*rho*T;


end
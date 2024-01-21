function [s1,t1,d1,x1,y1] = RefractRay(K,i,j,u,v)
if isfield(K,'UW')
    A = K.H*[i;j;u;v;1];
    s = A(1);
    t = A(2);
    x = A(3);
    y = A(4);
    nw = K.UW(1);
    nx = K.UW(2);
    ny = K.UW(3);
    nz = K.UW(4);
    nd = K.UW(5);
else
    s = K.ki*i;
    t = K.kj*j;
    x = K.ku*u + K.u0;
    y = K.kv*v + K.v0;
    nd = K.nd; nw = K.nw;
    nx = K.nx; ny = K.ny; nz = K.nz;
end
d0 = -nd/nz;
d1 = (d0*nz - s*nx - t*ny)/(nx*x + ny*y + nz);
s1 = s + d1*x;
t1 = t + d1*y;

dir0 = [x;y;1];
dir0N = sqrt(sum(dir0.^2));
dir0 = dir0./dir0N;
nVec = [nx;ny;nz];
NN = nx^2+ny^2+nz^2;
nVec = nVec./sqrt(NN);
dir1 = (1/nw)*cross(nVec,cross(-nVec,dir0)) - nVec.*sqrt(1-(1/nw)^2*sum(cross(nVec,dir0).*cross(nVec,dir0)));
dir1 = dir1./dir1(3);
x1 = dir1(1);
y1 = dir1(2);

end
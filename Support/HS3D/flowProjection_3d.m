function [u0,v0,w0] = flowProjection_3d( u,v,w, dx, dy, dz, Rp, Rpt )

     div = dx'*u + dy'*v+ dz'*w;     
     P = Rp\(Rpt\div);
     u0 = u - dx*P;
     v0 = v - dy*P;
     w0 = w - dz*P;
     
     
end
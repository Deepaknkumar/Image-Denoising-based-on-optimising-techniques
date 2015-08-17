function mxr = thresholdfuncsc4(u,xr)
    mxr = [];
            q = size(xr);
            thr = u(1);
            k = u(2);
            d = u(3);
            n = u(4);
             for m=1:q(2)
                if xr(1,m)>thr
                    mxr(1,m) = xr(1,m) -(((.5)*(thr^d)*(k))/(xr(1,m)^(d-1))) + (k-1)*thr;
                elseif abs(xr(1,m)) <= thr
                    mxr(1,m) = ((0.5)*(k*((abs(xr(1,m)))^n))*sign(xr(1,m)))/(thr^(n-1));
                elseif xr(1,m) < -thr
                    mxr(1,m) =xr(1,m) + (((.5)*((-thr)^d)*k)/(xr(1,m)^(d-1))) - (k-1)*thr;
                end
             end
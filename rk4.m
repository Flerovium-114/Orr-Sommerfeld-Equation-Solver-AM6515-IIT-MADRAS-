function [ t, w ] = rk4( fun, t0, tf, y0, h )
t = [t0:h:tf];
w = zeros(length(y0),length(t)); % initialize w array
w(:,1) = y0;
for i = 1:(length(t)-1);
    ti = t(i);
    wi = w(:,i);
    k1 = h*fun(ti,wi);
    k2 = h*fun(ti+h/2,wi+k1/2);
    k3 = h*fun(ti+h/2,wi+k2/2);
    k4 = h*fun(ti+h,wi+k3);
    w(:,i+1) = wi + 1/6*(k1+2*k2+2*k3+k4);
end

end


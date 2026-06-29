function [r,v,eunit]=eletorv(a,e,i,sig,w,f,myu)
p=a*(1-e^2);
h=sqrt(p*myu);
hunit=[sin(sig)*sin(i);-cos(sig)*sin(i);cos(i)];
eunit=[cos(w)*cos(sig)-cos(i)*sin(w)*sin(sig);cos(w)*sin(sig)+cos(i)*sin(w)*cos(sig);sin(w)*sin(i)];
ep=cross(hunit,eunit);
r=a*(1-e^2)/(1+e*cos(f))*(cos(f)*eunit+sin(f)*ep);
v=sqrt(myu/(a-a*e^2))*(-sin(f)*eunit+ep*(e+cos(f)));
end
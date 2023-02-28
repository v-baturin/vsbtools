function coef0=coef(f,v)
sizev=size(v);
sizef=size(f);
jj=1;
while jj<=max(sizef)
    ii=1;
    while ii<=max(sizev)
        coef1(jj,ii)=diff(sym(f(jj)),sym(v(ii)));
        ii=ii+1;
    end
    jj=jj+1;
    end
    coef0=eval(coef1);
end





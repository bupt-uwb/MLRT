function d=FrequencyCal(x,nw,ni)
n=nw;                                        
h=ni;                                       
s0=length(x);
win=hamming(n)';                            
c=1;
ncols=1+fix((s0-n)/h);                      
d=zeros((1+n/2),ncols);
for b=0:h:(s0-n)
    u=win.*x((b+1):(b+n));
    t=fft(u);
    d(:,c)=t(1:(1+n/2))';
    c=c+1;
end
end
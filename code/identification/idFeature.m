function [maxfHz,tcycle,tup,tdown,vup,vdown,upindensity,downindensity,ctwist, ...
    meanS,varS,skewnessS,kurtosisS,rmsS,fgfzS,arvS,pS,kaS,kfS] = idFeature(signal)
    fre=abs(fftshift(fft(signal,2^16)))*1000;
%     Y = (-2^16/2:2^16/2-1)*(20)/2^16;
%     figure(3)
%     plot(Y,fre);
%     xlabel('ÆµÂÊ/ºÕ×È');
%     ylabel('·ù¶È');
%     figure(4)
%     plot(signal)
    maxfA = max(fre(round(19/20*2^15):round(49/50*2^15)));
    maxfHz=10-((find(fre(round(19/20*2^15):round(49/50*2^15))==maxfA)+19/20*2^15)/(2^15)*10);
    c=size(signal);
    tcycle=c(1)-1;tup=0;tdown=0;
    twistStartT=1;twistStartA=signal(1);upOrdown=0;Lvup=[];Lvdown=[];
    for i = 2:tcycle
        if(signal(i)>signal(i-1))
            tup = tup+1;
        else
            tdown = tdown+1;
        end
        if(i==2)
            if(signal(i)>signal(i-1))
                upOrdown=1;
            else
                upOrdown=2;
            end
        else
            if(signal(i)>signal(i-1) && upOrdown==1)
                continue
            elseif(signal(i)>signal(i-1) && upOrdown==2)
                upOrdown=1;
                Lvdown = [Lvdown,(signal(i-1)-twistStartA)/(i-1-twistStartT)];
                twistStartT=i-1;twistStartA=signal(i-1);
            elseif(signal(i)<signal(i-1) && upOrdown==2)
                continue
            elseif(signal(i)<signal(i-1) && upOrdown==1)
                upOrdown=2;
                Lvup = [Lvup,(signal(i-1)-twistStartA)/(i-1-twistStartT)];
                twistStartT=i-1;twistStartA=signal(i-1);
            else
                fprintf('warning!!!!!!!!!!!!!!!!!!!!!!!!!');
            end
        end
    end
    vup = mean(Lvup);vdown = mean(Lvdown);
    upindensity = max(signal);downindensity=min(signal);
    upsize =size(Lvup);downsize =size(Lvdown);
    ctwist = min(upsize(2),downsize(2));
    
    meanS = mean(signal); varS = var(signal);skewnessS=skewness(signal);kurtosisS = kurtosis(signal);
    rmsS = rms(signal);fgfzS = mean(sqrt(abs(signal))).^2;arvS = mean(abs(signal));
    pS = max(signal) - min(signal);kaS = rmsS./arvS; kfS = pS./fgfzS;
    
%     meanS = mean(signal); varS = var(signal);skewnessS=skewness(signal);kurtosisS = kurtosis(signal);
%     rmsS = rms(signal);fgfzS = mean(sqrt(abs(signal)),2).^2;arvS = mean(abs(signal));
%     pS = max(signal) - min(Signal);kaS = rmsS/arvS; kfS = pS/fgfzS;
    %maxfHz,tcycle,tup,tdown,vup,vdown,upindensity,downindensity£¬ctwist
%
end
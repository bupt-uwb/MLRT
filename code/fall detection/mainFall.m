clc
close all
clear
path = 'H:\data\20210511\20210607\';
person = 'cw';
for number = 3:30
    close all
    version = 'Id8';
    radar1 = load([path,person,'_fall',num2str(number),'.mat']);
    % radar1 = load('H:\data\20210511\Data2\20210511_id_dy_sit.mat');
    radar1=radar1.data;
    radar1(:,594)=[];
    radar1(:,593)=[];
    radar1(:,592)=[];
    RawData = radar1;

    % PureData = RawData-mean(RawData);
    PureData = newfilter(radar1(:,:));
    % PureData = PureData(:,50:end);
%     figure(1)
%     mesh(abs(PureData))
%     title('PureData')
    c=size(PureData);
    for i = 1:c(1)
        maxV(i) = max(abs(PureData(i,:)));
        lmaxRange(i)=(find(abs(PureData(i,:))==maxV(i)));
    end
    for ii = 1:c(1)-20
        leiji(ii) = sum(sum(abs(PureData(ii:ii+19,200:end))));
    end

    picS = sum(abs(PureData(:,200:end)),2);
    xtk_label = {'1s','2s','3s','4s','5s','6s','7s','8s','9s','10s'};
%     ytk_label={'0',''};
    figure(1)
    plot(picS)
    title('The time series')
    
    
    set(gca,'xticklabel',xtk_label);
%     set(gca,'yticklabel',ytk_label);
    xlabel('Time');
    ylabel('Amplitude');
    set(gca, 'Fontname','times new roman', 'Fontsize', 20,'FontWeight','bold')
    
    N=380;  
    Xkf=zeros(1,N);
    Xkf(1)=lmaxRange(1);
    F=1;
    G=1;
    H=1;
    I=1; 
    P(1)=0.01;
    Q=0.01;
    R=1;
    %%%%%%%%%%%%%%%%%%%%%%%
    for k=2:N 
        X_pre=F*Xkf(k-1);           
        P_pre=F*P(k-1)*F'+Q;        
        Kg=P_pre*inv(H*P_pre*H'+R); 
        e=lmaxRange(k)-H*X_pre;            
        Xkf(k)=X_pre+Kg*e;         
        P(k)=(I-Kg*H)*P_pre;
    end
    t=1:N;
    figure(2);
    Mrange= ((lmaxRange+50)*0.642978+20)/100;
    MXkf = ((Xkf+50)*0.642978+20)/100;
    plot(t,Mrange,'r',t,MXkf,'b')
    legend('measure','kalman extimate');
    title('Trajectory of the target')
    Range_label = {'0','0.5m','1m','1.5m','2m','2.5m','3m','3.5m','4m'};
    axis([0,180,0,4])
    set(gca,'xticklabel',xtk_label);
    set(gca,'yticklabel',Range_label);
    xlabel('Time');
    ylabel('Range');
    
%     leiji = ones(1,380);
%     for ii = 1:c(1)
%             buffer = abs(PureData(ii,lmaxRange(ii)));
%             leiji(ii) = abs(PureData(ii,lmaxRange(ii)));
%     end
    
%     figure(2)
%     plot(maxV)
%     title('maxV')
%     figure(3)
%     plot(leiji)
%     title('leiji')
%     cSignal = smoothdata(maxV,'gaussian',5);
%     figure(3)
%     plot(cSignal)
%     title('cSignal')
%     fSignal = filter(ones(1, 5)/5, 1, maxV);
%     figure(3)
%     plot(fSignal)
%     title('fSignal')
    for j =1:c(1)-40
        pe(j,:) = pec(maxV(1,j:j+19),3,3);
    %     se(j,:) = multiscaleSampleEntropy_compatible(maxV(1,j:j+19));%0.1~0.25
        ae(j,:) = ApEn(3,4e-2,leiji(1,j:j+19),2);
%         ae(j,:) = ApEn(3,6e-4,leiji(1,j:j+19),2);
    %     pe(j,:) = MPerm(maxV(1,j:j+19),3,1,2);
    end
    for j =1:c(1)-20
       
    %     se(j,:) = multiscaleSampleEntropy_compatible(maxV(1,j:j+19));%0.1~0.25
        aee(j,:) = ApEn(3,1e-1,leiji(1,1:j),5);
    %     pe(j,:) = MPerm(maxV(1,j:j+19),3,1,2);
    end
%     figure(4)
%     plot(pe)
%     title('pe')
    figure
    xtk_label = {'2s','3s','4s','5s','6s','7s','8s','9s'};
    plot(2:0.05:18.95,ae)
    title('The set of ApEn in time series')
    set(gca,'xticklabel',xtk_label);
%     set(gca,'yticklabel',ytk_label);
    xlabel('Time');
    ylabel('Value');
    
%     figure(6)
%     plot(aee)
%     title('aee')
    if(min(ae(201:end))<-0.13)%min(pe(70:end))<1.2
        
        fprintf(['fall ',num2str(number),...
            ' ',num2str(min(ae(201:end))),' ',num2str(min(find(ae(201:end)==(min(ae(201:end)))))/20+12.4),'\n'])
    else
        fprintf(['non-fall ',num2str(number),' ',num2str(min(pe(201:end))),...
            ' ',num2str(min(ae(201:end))),'\n'])
    end
    
end
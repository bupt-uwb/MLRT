clc
close all
clear
path = 'H:\data\20210511\20210607\';
person = 'cwh';
version = 'id8';
confeature=[];
for number = 27:30
%     close all
    number
%     close all
    radar1 = load([path,person,'_fall',num2str(number),'.mat']);
    % radar1 = load('H:\data\20210511\Data2\20210511_id_dy_sit.mat');
    radar1=radar1.data;
    radar1(:,594)=[];
    radar1(:,593)=[];
    radar1(:,592)=[];

    % PureData = RawData-mean(RawData);
    PureData = newfilter(radar1(:,:));
    % PureData = PureData(500:5000,:);
    % mesh(PureData)
    % title('PureData')
%     waveform = phased.LinearFMWaveform('PulseWidth',1e-3,'PRF',23.328e7/2,'OutputFormat','Pulses','NumPulses',1,'SweepBandwidth',1e5);
%     prf = waveform.PRF;
%     fs =prf*200;
%     response = phased.RangeDopplerResponse('DopplerFFTLengthSource','Property','DopplerFFTLength',2048,'SampleRate',fs,'DopplerOutput','Speed','OperatingFrequency',7.29e9,'PRFSource','Property','PRF',prf);
%     filt = getMatchedFilter(waveform);
%     [resp,rng_grid,dop_grid] = response(PureData(1:200,:)',filt);
%     figure
%     plotResponse(response,PureData(1:200,:)',filt,'Unit','db')
% %     ylim([0 12000])         
% 
%     fasttime = unigrid(0,1/fs,1/prf,'[)');
%     rangebins = (physconst('Lightspeed')*fasttime/2);
%     figure
%     plot(rangebins,abs(PureData(1:200,:)'))

    cPure = size(PureData);
    for i = 1:200
        maxV(i) = max(abs(PureData(i,:)));
        lmaxRange(i)=(find(abs(PureData(i,:))==maxV(i)));
    end

    N=200;  
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
    figure(1);
    Mrange= ((lmaxRange+50)*0.642978+20)/100;
    MXkf = ((Xkf+50)*0.642978+20)/100;
    plot(t,Mrange,'r',t,MXkf,'b')
    legend('measure','kalman extimate');
    title('Trajectory of the target')
    xtk_label = {'1s','2s','3s','4s','5s','6s','7s','8s','9s','10s'};
    Range_label = {'0','0.5m','1m','1.5m','2m','2.5m','3m','3.5m','4m'};
    axis([0,180,0,4])
    set(gca,'xticklabel',xtk_label);
    set(gca,'yticklabel',Range_label);
    xlabel('Time');
    ylabel('Range');
    for i2 = 1:200
        trackData(i2) = PureData(i2,round(Xkf(1,i2)));
    end
    
    figure(2);
    plot(trackData);
    title('trackData')
    move50Data = ones(200,50);
    
    for i1  = 1:200
        if(round(Xkf(1,i1))+25 > 541)
            move50Data(i1,:) = PureData(i1,541-49:541);
        else
            move50Data(i1,:) = PureData(i1,round(Xkf(1,i1))-24:round(Xkf(1,i1))+25);
%         round(Xkf(1,i1))
        end
    end
    
    figure(3)
    mesh(PureData)
    title('PureData')
    testsignal = move50Data(:,25);

    testsignalinterval = move50Data;

%     psd = abs(fftshift(fft2(testsignalinterval'))).^2;
%     % 通过对数变换，便于观察
%     psd = 10 * log10(psd); 
%     figure
%     imagesc(psd)
%     mesh(psd)

    figure(4);
    fre=abs(fftshift(fft(testsignal,2^16)))*1000;
    Y = (-2^16/2:2^16/2-1)*(20)/2^16;
    plot(Y,fre);
    title('F testsignal')
    xlabel('频率/赫兹');
    ylabel('幅度');

    % 使用二项分布的系数作为权值
    % 加权移动平均滤波器
    h = [1/2 1/2];
    binomialCoeff = conv(h,h);
    for n = 1:4
        binomialCoeff = conv(binomialCoeff,h);
    end
    binomialMA = filter(binomialCoeff, 1, testsignal);
    %对数滤波
    alpha = 0.50;
    exponentialMA = filter(alpha, [1 alpha-1], testsignal);
    %平滑
    c = smoothdata(testsignal,'gaussian',20);
    
    fre=abs(fftshift(fft(c,2^16)))*1000;
    figure(5)
    plot(Y,fre);
    title('F c')
    xlabel('频率/赫兹');
    ylabel('幅度');

    b=findMin(c);
    % b = 1:100:size(c);
    figure(6);hold on;box on;
    title('c and b')
    % plot(1:size(binomialMA),binomialMA,'black');
    plot(1:N,c,'b',b,c(b),'r^');
    % plot(b,binomialMA(b),'r')
    xtk_label={'100s','105s','110s','115s','120s','125s','130s','135s','140s','145s','150s'};
    % axis([2000,3000,-4e-3,4e-3])
    legend('The filtered signal','Break point')%'Hourly Temp', ...'Exponential Weighted Average','location','best'
    ylabel('Amplitude')
    xlabel('Time')
    set(gca,'xticklabel',xtk_label);
    title('Signal filtered by Gaussian window')

    set(gca,'xticklabel',xtk_label);

    nfeature=[];
    nfeatureInterval=[];
    % 
    cc = size(b);
    if(cc(1)==1)
        csize = cc(2);
    else
        csize = cc(1);
    end
    feature=[];
%     csize
    for ir = 1:180
        rangef(ir) = sum(MXkf(ir:ir+19))/20;
    end
    Vel = diff(rangef,1)*20;
    acce = diff(rangef,2)*20;
    figure(7)
    plot(rangef)
    title('range')
    figure(8)
    plot(Vel)
    title('Vel')
    figure(9)
    plot(acce)
    title('acce')
    figure(10)
    plot(rangef(1:179),Vel,'o')
    title('range-v')
    figure(11)
    plot(rangef(1:178),acce,'o')
    title('Distance-Acceleration','FontWeight','bold')
    ylabel('Acceleration(m/s^2)','FontWeight','bold')
    xlabel('Distance(m)','FontWeight','bold')
    axis([2.7,3.3,-0.05,0.05])
    set(gca,'FontSize',20)
    set(gca,'XTick',2.7:0.1:3.3);
    set(gca,'YTick',-0.05:0.02:0.05);
    
%     xlabel('Frequency (Hz)','FontWeight','bold')
%     ylabel('Amplitude','FontWeight','bold')
    set(gca, 'Fontname','times new roman', 'Fontsize', 20,'FontWeight','bold')
%     set(gca,'FontSize',20)

    figure(12)
    plot(Vel(1:178),acce,'o')
    title('Velocity-Acceleration','FontWeight','bold')
    ylabel('Acceleration(m/s^2)','FontWeight','bold')
    xlabel('Velocity(m/s)')
    axis([-0.4,0.4,-0.05,0.05])
    set(gca,'FontSize',20)
    set(gca,'XTick',-0.4:0.1:0.4);
    set(gca,'YTick',-0.05:0.02:0.05);
    set(gca, 'Fontname','times new roman', 'Fontsize', 20,'FontWeight','bold')
%     [signal,fsc] = wavread('tone4.wav');
%     nw=16;ni=nw/4;
%     d=FrequencyCal(testsignal',nw,ni);
%     tt=(0:ni:(length(testsignal)-nw))/10;
%     ff=(0:(nw/2))*10/nw*2;
%     % imagesc(tt,ff,20*log10(abs(d)));
%     figure
%     imagesc(tt,ff,abs(d).^2);
%     xlabel('Time(s)');
%     ylabel('Frequency(Hz)')
%     title('时频图')
%     axis xy
    rangeAcce = [rangef(1:178);acce];
    VAcce = [Vel(1:178);acce];
    RVA = [rangef(1:178);Vel(1:178);acce];
%     X1 = rangeAcce - repmat(mean(rangeAcce),2,1);
%     X2 = VAcce - repmat(mean(VAcce),2,1);
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED]=pca(RVA');
%     COEFF(:,1)
%     LATENT
    var(acce)
    eva1 = evalclusters(rangeAcce','kmeans','CalinskiHarabasz','KList',1:10);
    eva2 = evalclusters(VAcce','kmeans','CalinskiHarabasz','KList',1:10);
%     figure(13)
%     [idx,C] = kmeans(rangeAcce',2);
%     idx2Region = kmeans(rangeAcce',2,'Start',C);
%     gscatter(rangef(1:178),acce,idx2Region)
%     legend('Region 1','Region 2');
%     figure(14)
%     [idx,C] = kmeans(VAcce',2);
%     idx2Region = kmeans(VAcce',2,'Start',C);
%     gscatter(Vel(1:178),acce,idx2Region)
%     legend('Region 1','Region 2');

    movefeature = [COEFF(:,1);LATENT;mean(RVA(1,:));var(RVA(1,:));mean(RVA(2,:));var(RVA(2,:));mean(RVA(3,:));var(RVA(3,:))];
    
    %     opts = statset('Display','final');
%     [idx,ctrs] = kmeans(rangeAcce',2,...
%     'Distance','city',...
%     'Replicates',5);
%     figure
% %     plot(rangeAcce(idx==1,1),rangeAcce(idx==1,2),'r.','MarkerSize',12)
% %     hold on
% %     plot(rangeAcce(idx==2,1),rangeAcce(idx==2,2),'m.','MarkerSize',12)
%     
%     plot(ctrs(1,:),ctrs(2,:),'ko',...
%     'MarkerSize',2,'LineWidth',1.5)
%     legend('Cluster 1','Cluster 2','Centroids',...
%     'Location','NW')


    if  csize>1
        for i = 2:csize
         [maxfHz,tcycle,tup,tdown,vup,vdown,upindensity,downindensity,ctwist, ...
        meanS,varS,skewnessS,kurtosisS,rmsS,fgfzS,arvS,pS,kaS,kfS] = idFeature(c(b(i-1):b(i)));
        
        oneFeature =  [maxfHz,tcycle,tup,tdown,vup,vdown,upindensity,downindensity,ctwist, ...
        meanS,varS,skewnessS,kurtosisS,rmsS,fgfzS,arvS,pS,kaS,kfS];
        [DDmeanS,DDvarS,DDskewnessS,DDkurtosisS,DDrmsS,DDfgfzS,DDarvS,DDpS,DDkaS,DDkfS,...
        DVmeanS,DVvarS,DVskewnessS,DVkurtosisS,DVrmsS,DVfgfzS,DVarvS,DVpS,DVkaS,DVkfS] = idFeatureInterval(testsignalinterval(b(i-1):b(i),:));
        
        oneFeatureInterval=[DDmeanS,DDvarS,DDskewnessS,DDkurtosisS,DDrmsS,DDfgfzS,DDarvS,DDpS,DDkaS,DDkfS,...
        DVmeanS,DVvarS,DVskewnessS,DVkurtosisS,DVrmsS,DVfgfzS,DVarvS,DVpS,DVkaS,DVkfS];
        nfeature = [nfeature;oneFeature];
        nfeatureInterval = [nfeatureInterval;oneFeatureInterval];
        
        
        nsumfeature = [nfeature,nfeatureInterval]; 
        nsumfeature = mean(nsumfeature,1);
        end
    else
        fprintf('!!!!!!!!!!!!!!!!!!!!')
    end    
%     size(nsumfeature)
%     feature = [nfeature,nfeatureInterval];
    if isempty(confeature)
        confeature = [nsumfeature';movefeature];        
    else 
        confeature = [confeature,[nsumfeature';movefeature]];
    end
    
end
savename = ['H:\data\20210511\feature\',version,'\',person,'FId','.mat'];
save(savename,'confeature')
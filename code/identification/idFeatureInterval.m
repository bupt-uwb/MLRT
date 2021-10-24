function [DDmeanS,DDvarS,DDskewnessS,DDkurtosisS,DDrmsS,DDfgfzS,DDarvS,DDpS,DDkaS,DDkfS,...
    DVmeanS,DVvarS,DVskewnessS,DVkurtosisS,DVrmsS,DVfgfzS,DVarvS,DVpS,DVkaS,DVkfS] = idFeatureInterval(signal)
    meanS = mean(signal); varS = var(signal);skewnessS=skewness(signal);kurtosisS = kurtosis(signal);
    rmsS = rms(signal);fgfzS = mean(sqrt(abs(signal))).^2;arvS = mean(abs(signal));
    pS = max(signal) - min(signal);kaS = rmsS./arvS; kfS = pS./fgfzS;
    
    DDmeanS = mean(meanS);DDvarS = mean(varS);DDskewnessS = mean(skewnessS);DDkurtosisS = mean(kurtosisS);
    DDrmsS = mean(rmsS);DDfgfzS = mean(fgfzS);DDarvS = mean(arvS);DDpS = mean(pS);
    DDkaS = mean(kaS);DDkfS = mean(kfS);
    
    DVmeanS = var(meanS);DVvarS = var(varS);DVskewnessS = var(skewnessS);DVkurtosisS = var(kurtosisS);
    DVrmsS = var(rmsS);DVfgfzS = var(fgfzS);DVarvS = var(arvS);DVpS = var(pS);
    DVkaS = var(kaS);DVkfS = var(kfS);
end

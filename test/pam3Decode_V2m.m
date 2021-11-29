% cut_off_frequency = 125e6;
% verbose = 1;
% if(verbose); tic;end
% [X,Y] = interpolateTo1GHzAndFilter(scopeObject,6,5,cut_off_frequency,0);
% if(verbose);disp(['Interpolating to 1 GHz and filtering execution time: ' num2str(toc) newline]);end 
% % X = scopeObject.time;
% % Y= scopeObject.channels(6).value - scopeObject.channels(5).value;
% Ytest = Y/max(Y);
[Xt,Yt] = pam3(scopeObject,X,Ytest,verbose);
PNchannels =2:3;
[Ra, Rb] = demultiplexer(Yt);

i=1;
j=1;
successiveIndices=0;
ssdIndex = 1;
prevIndex = 0;
threeBit = [];
ssdIndices = [];
esdIndices = [];
index=0;
while index <numel(Ra)
    index = index +1;
    threeBit = [threeBit; twoTernaryToDataSymbol(Ra(index),Rb(index))];
      
    checkSSD = isnan(threeBit);
    if  numel(threeBit)==9 
        if sum(checkSSD)==9 % SSD code group {00},{00},{00} with 00 replaced by nan,nan,nan
            ssdIndices(i) = index*2;
            i=i+1;
            threeBit = [];
        elseif sum(checkSSD)==6 && sum(threeBit(end-2:end))==3 % ESD code group {00},{00},{11}
            esdIndices(j) = index*2;
            j=j+1;
            threeBit = [];
%         elseif sum(checkSSD)==3 && sum(threeBit(4:6))==3 && sum(threeBit(1:3))== 0 % shifted symbol {00},{11},{-1-1}
%             disp('Symbols shifted');
%             Ra(2:end)=Ra(1:end-1);
%             index =1;
%             j=j+1;
%             threeBit = [];  
        end
    end
    % Remove the first three bits to keep 9 bits
    if numel(threeBit)>=9
        threeBit(1:3)=[];
    end
    
    prevIndex = index;
end

figure;
for i= 1 : numel(PNchannels)
    scopeObject.channels(PNchannels(i)).pn.plot(0,10-(i-1)*3,1,'g');
end
hold on
for i=1 : numel(ssdIndices(1:100))
    xline(scopeObject.time(Xt(ssdIndices(i))),'c--');
end

figure
startIndex = 56200;
endIndex = 56500;
plot(X(Xt(startIndex):Xt(endIndex)),Y(Xt(startIndex):Xt(endIndex)))
hold on
plot(X(Xt(startIndex:endIndex)),Yt(startIndex:endIndex)/10,'*')
xline(scopeObject.channels(2).pn(1, 1).time,'g--')
xline(scopeObject.channels(3).pn(1, 1).time,'g--')

for i=Xt(startIndex:endIndex)
    xline(X(i));
end

yline(0.05393,'r')
yline(-0.057792,'r')
ylim([-1.1 1.1])
% hold on
% for i=1 : numel(esdIndices)
%     xline(scopeObject.time(Xt(esdIndices(i))),'r--');
% end
     
% findSSD = find(Yt ==[0,0]);
% indexTa = 1:2:numel(Yt);
% indexTb = 2:2:numel(Yt);
% Ta = Yt(indexTa);
% Tb = Yt(indexTb);
% findSSD = ismember(Yt,[0 0]);
% locationSSD = find(findSSD==1);
% Yt(locationSSD(1)-5:locationSSD(1)+5)
function [objScope, chNr, threshold, cut_off_frequency, verbose] = splitVarargin(varargin)
    varargin = varargin{1};
    if(numel(varargin) >= 1)
        if(isa(varargin{1}, 'scope'))
            objScope = varargin{1};
        else
            ethPHYdecode.help;
            error('Need a scope object as first input');
        end
        if(numel(varargin) >= 2)
            if(isa(varargin{2}, 'double') && numel(varargin)==2)
                verbose = varargin{2};
            else
                varargin(1) = [];
                while ~isempty(varargin)
                    if(isstruct(varargin{1}))
                        var = varargin{1};
                        fields = fieldnames(var);
                        varargin(length(fields)*2+1:end+length(fields)*2-1)=varargin(2:end);
                        for i = 1:numel(fields)
                            varargin(i*2-1)=fields(i);
                            varargin{i*2}=var.(char(fields(1)));
                        end
                    elseif(ischar(varargin{1}))
                        switch lower(varargin{1})
                            case 'verbose'
                                verbose = varargin{2};
                                varargin(1:2) = [];
                            case 'threshold'
                                threshold = varargin{2};
                                varargin(1:2) = [];
                            case 'channelnr'
                                chNr = varargin{2};
                                varargin(1:2) = [];
                            case 'cut_off_frequency'
                                cut_off_frequency = varargin{2};
                                varargin(1:2) = [];
                            otherwise
                                warn('Unknown argument');
                                varargin(1) = [];
                        end
                    else
                        warn('Unknown argument');
                        varargin(1) = [];
                    end
                end
            end
        end
    end

    if(~exist('verbose','var'));verbose=-1;warn('All underlying functions are executed in verbose mode');end
    if(~exist('threshold','var'));threshold=0.5;end
    if(~exist('cut_off_frequency','var'));cut_off_frequency=2*125e6;end
end

function [adjustedTime,adjustedSignal] = interpolateTo1GHzAndFilter(objScope,plusChannel,minChannel,cut_off_frequency,verbose)
scale=1;
freq = 8e-9;

X = objScope.time;
adjustedSignal = objScope.channels(plusChannel).value - objScope.channels(minChannel).value;

while objScope.sample_interval > 1e-9
    objScope.sample_length = objScope.sample_length * 2;
    objScope.sample_interval = objScope.sample_interval / 2;
    scale = scale * 2;
end

adjustedTime = X(1)-objScope.sample_interval*(scale-3)/2:objScope.sample_interval:X(end)+objScope.sample_interval*(scale+1)/2;
expfft = ceil(log2(length(adjustedSignal)));
F = fft(adjustedSignal,2^expfft);
Forig = F;

% Adding zero's (resulting in 1 GS/s)
F = [F(1:floor((end+1)/2))*scale, zeros(1,length(F)*(scale-1)), F(ceil((end+1)/2):end)*scale];
Fs = 1/objScope.sample_interval; % Sampling frequency
L = objScope.sample_length;      % Length of signal
% Cutting 40-60Hz, 140-160Hz and Cut-off frequencies
filterfreq = [40,60;140,160;240,260;340,360;cut_off_frequency,1e99];
for i = 1:size(filterfreq,1)
    locmin([1,2]) = round(filterfreq(i,1)*L/Fs)+1;
    locmax([1,2]) = round(filterfreq(i,2)*L/Fs)+1;
    if(locmin(1)<1);locmin(1)=1;end;
    if(locmax(1)<1);locmax(1)=1;end;
    if(locmin(2)<2);locmin(2)=2;end;
    if(locmax(2)<2);locmax(2)=2;end;
    if(locmin(1)>length(F)/2);locmin(1)=round(length(F)/2);end;
    if(locmax(1)>length(F)/2);locmax(1)=round(length(F)/2);end;
    if(locmin(2)>length(F)/2);locmin(2)=round(length(F)/2);end;
    if(locmax(2)>length(F)/2);locmax(2)=round(length(F)/2);end;
    L1 = locmin(1):locmax(1);
    L2 = length(F)-(locmax(2)-2):length(F)-(locmin(2)-2);
    F(L1) = 0;
    F(L2) = 0;
end
Yorig = adjustedSignal;
adjustedSignal = real(ifft(F,(2^expfft)*scale));
adjustedSignal = adjustedSignal(1:length(adjustedTime));

%% SETTING PLOT PARAMETERS
if(verbose)
    plotMin =  500;
    plotMax =  1000;
    Xmin = X(plotMin);
    Xmax = X(plotMax);
    Ymin = -1.5;
    Ymax = +1.5;
end
if(verbose)
    figure;
    hold on;
    Xorig = objScope.time;
    plt(Xorig(Xorig>Xmin & Xorig<Xmax),Yorig(Xorig>Xmin & Xorig<Xmax),'downsample',1e5,'color', [0.25 0.25 0.25]);
    plt(X(plotMin:plotMax),adjustedSignal(plotMin:plotMax),'downsample',1e5,'r');
    xlabel('t(s)');
    ylabel('U(V)');
    axis([Xmin,Xmax,Ymin,Ymax]);
    hold off;
    
    figure;
    f = Fs*(0:ceil(L/2))/L;
    P1 = abs(F(1:ceil(L/2)+1));
    %ax = [1 max(f) min(P1(P1>10^(-10))) 1];
    plt(f(f<256e6)/2/1e6,P1(f<256e6),color.ch2);
    axis([0 125 0  15e5]);
    
    xlabel('f [MHz])');
    
    %                 axis(ax);
    %                 subplot(2,1,2);
    %                 P2 = P2(1:ceil(L/2)+1);
    %                 loglog(f,P1,color.lightgrey);
    %                 axis(ax);
    %                 hold on;
    %                 loglog(f,P2,color.ch2);
    %                 title('Single-Sided Amplitude Spectrum of P2(t)');
    %                 xlabel('f (Hz)');
    %                 ylabel('|P2(f)|');
    %                 axis(ax);
    %                 hold off;
end
%% PLOTTING GRAPH 1
if(verbose)
    figure;
    subplot(4,1,1);
    XT = objScope.time;
    YT = objScope.channels(1).value - objScope.channels(2).value;
    
    hold on;
    plot(XT(Xmin<XT & XT<Xmax),YT(Xmin<XT & XT<Xmax), 'color', [0.5 0.5 0.5]);
    plot(adjustedTime(plotMin:plotMax),adjustedSignal(plotMin:plotMax),'g');
    axis([Xmin,Xmax,Ymin,Ymax]);
    hold off;
end
% objScope.time = X;
clear F Forig Yorig Scale filterfreq;
clear X expfft
end

function[ternarySampleHitIndex, ternarySample]= pam3(objScope,timeScope, signal,verbose)
    if (verbose);tic; end
    % 100BASE-T1
    Ft = 66.66e6;
    ThresholdCondition = false;
    percentOnes = 0;
    percentMinusOnes = 0;
    posThresholdOffset = 0;
    negThresholdOffset = 0;
    % calculate threshold indices
    posSignal = nan(1, numel(signal));
    posSignalIndex = signal >= max(signal)*0.2 & signal <= (max(signal)*0.6);
    posSignal(posSignalIndex)= signal(posSignalIndex);
    posPeaks = findpeaks(posSignal);
    
    negSignal = nan(1, numel(signal));
    negSignalIndex = signal<min(signal)*0.1 & signal>=(min(signal)*0.8);
    negSignal(negSignalIndex)= signal(negSignalIndex);
    negPeaks = findpeaks(negSignal);
    
    
        Ytest = Y/max(Y);
    Ydiff = [diff(Ytest),0];
    Xtest = repelem(Xt,5);
    Xtest(1:5:end)= Xt-2;
    Xtest(2:5:end)= Xt-1;
    Xtest(4:5:end)= Xt+1;
    Xtest(5:5:end)= Xt+2;
    
    
    while ~ThresholdCondition      
        % Calculate  thresholds
        posThreshold = min(posPeaks) + posThresholdOffset; 
        negThreshold = mean(negPeaks)+ negThresholdOffset;
  
        % Search indices of plus and minus ones
        plusOneIndex = signal >= posThreshold;
        minusOneIndex = signal <= negThreshold;

        % Find zero indices
%         diffSignal = [0, diff(signal)];
%         DiffThresholdZero = max(abs(diffSignal))/2; 
        zeroIndex = signal < posThreshold & signal > negThreshold ; %& (diffSignal < DiffThresholdZero & diffSignal > -DiffThresholdZero)

        % decode to ternary bitstream
        ternaryBit = nan(1,length(signal));
        ternaryBit(plusOneIndex) = 1;
        ternaryBit(minusOneIndex) = -1;
        ternaryBit(zeroIndex) = 0;

        % Calculate Ternary bit on sample hits
        firstPeakID = find(signal==findpeaks(signal, 'NPeaks',1),1,'first');

        ternarySampleHitIndex = timeScope(firstPeakID):(1/Ft):timeScope(end);
        ternarySampleHitIndex=round(ternarySampleHitIndex./objScope.sample_interval)-1;
        ternarySample = ternaryBit(ternarySampleHitIndex);

        % calculate distribution
        symbolLength = numel(ternarySample);
        amountOfOnes = sum(ternarySample==1);
        amountOfZeros = sum(ternarySample==0);
        amountOfMinusOnes = sum(ternarySample==-1);
        amountOfNaNs = sum(isnan(ternarySample));
        percentOnes = amountOfOnes/symbolLength*100;
        percentZeros = amountOfZeros/symbolLength*100;
        percentMinusOnes = amountOfMinusOnes/symbolLength*100;  
        
        if verbose           
            disp('----- Threshold levels -----')
            disp(['Positive threshold: '  num2str(posThreshold)])
            disp(['negative threshold: ' num2str(negThreshold) newline])
            
            disp('------ PAM3 information ------');        
            disp(['# ones: '  num2str(amountOfOnes) ' (' num2str(round(percentOnes,2)) '%)']);
            disp(['# zeros: '  num2str(amountOfZeros) ' (' num2str(round(percentZeros,2)) '%)']);
            disp(['# minus ones: ' num2str(amountOfMinusOnes) ' (' num2str(round(percentMinusOnes,2)) '%)']);
            disp(['# NaN: ' num2str(amountOfNaNs) ' (' num2str(round(amountOfNaNs/symbolLength*100,2)) '%)' newline]);
        end
        
        % check percentages are between 25 and 35 %
        if  checkPercentage(percentOnes)==0 && checkPercentage(percentMinusOnes)==0 
            ThresholdCondition = true;
        else % if not adjust positive and negative threshold independently         
            % check positive threshold
            if  checkPercentage(percentOnes)==1
                posThresholdOffset =  posThresholdOffset + min(posPeaks)*0.1;
            elseif checkPercentage(percentOnes)==-1
                posThresholdOffset =  posThresholdOffset - min(posPeaks)*0.1;
            end
            
            % check negative threshold
            if  checkPercentage(percentMinusOnes)==1
                negThresholdOffset =  negThresholdOffset  - abs(min(negPeaks))*0.1;
            elseif checkPercentage(percentMinusOnes)==-1
                negThresholdOffset =  negThresholdOffset + abs(min(negPeaks))*0.1;
            end         
        end   
    end
    
    if verbose
        figure
        plot(timeScope(1:1000),ternaryBit(1:1000))
        hold on
        plot(timeScope(1:1000),signal(1:1000),'-o')
        yline(posThreshold,'--g','+Threshold')
        yline(negThreshold,'--g','-Threshold')

        figure
        hold on
        plot(timeScope(1:4000),signal(1:4000))
        ternaryPlotSize= sum(timeScope(ternarySampleHitIndex)<=timeScope(4000));
        plot(timeScope(ternarySampleHitIndex(1:ternaryPlotSize)),ternarySample(1:ternaryPlotSize)/10,'*');
%         plot(timeScope(1:4000),diffSignal(1:4000))
        legend('signal','ternary bit','diff signal')
    %     plot(timeScope(firstPeakID),signal(firstPeakID),'go')

        for i=1 : find(timeScope(ternarySampleHitIndex)<= timeScope(4000),1,'last')
            xline(timeScope(ternarySampleHitIndex(i)));
        end

        yline(posThreshold,'--g','+Threshold')
        yline(negThreshold,'--g','-Threshold')
   end
    
    % local function to check the percentage is between 25 and 35 %
    function check = checkPercentage(percentage)
        if percentage >35 
            check = 1;
        elseif  percentage <25 
            check = -1;
        else
            check = 0;
        end
    end

    clear symbolLength amountOfOnes amountOfZeros amountOfMinusOnes amountOfNaNs firstPeakID ternaryBit
    if(verbose);disp(['PAM3 execution time: ' num2str(toc) newline]);end      
end

function[Ra, Rb]= demultiplexer(ternarySample)
Ra(:) = int16(ternarySample(1:2:(numel(ternarySample)-1)));
Rb(:) = int16(ternarySample(2:2:numel(ternarySample)));
end

function ScrambleBits = twoTernaryToDataSymbol(Ra,Rb)
    switch true
        case Ra == -1 && Rb == -1 
            ScrambleBits = [0;0;0];
        case Ra == -1 && Rb == 0
            ScrambleBits = [0;0;1];
        case Ra == -1 && Rb == 1
            ScrambleBits = [0;1;0];
        case Ra == 0 && Rb == -1
            ScrambleBits = [0;1;1];
        case Ra == 0 && Rb == 0
            ScrambleBits = [nan;nan;nan];
        case Ra== 0 && Rb == 1
            ScrambleBits = [1;0;0];    
        case Ra == 1 && Rb == -1
            ScrambleBits = [1;0;1];
        case Ra == 1 && Rb == 0
            ScrambleBits = [1;1;0];   
        case Ra == 1 && Rb == 1
            ScrambleBits = [1;1;1];    
        otherwise
            disp('data symbol mapping failed');
    end 
end
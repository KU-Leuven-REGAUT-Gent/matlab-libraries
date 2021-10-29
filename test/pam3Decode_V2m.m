cut_off_frequency = 125e6;
verbose = -1;
[X,Y] = interpolateTo1GHzAndFilter(scopeObject,cut_off_frequency,0);
[Xt,Yt] = pam3(scopeObject,X,Y);



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

function [adjustedTime,adjustedSignal] = interpolateTo1GHzAndFilter(objScope,cut_off_frequency,verbose)
scale=1;
freq = 8e-9;

X = objScope.time;
adjustedSignal = objScope.channels(1).value - objScope.channels(2).value;

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

function[ternaryTime, ternaryBit]= pam3(objScope,time, signal)

    % 100BASE-T1
    Ft = 66.66e6;

    % calculate threshold indices
    threshold=max(abs(signal))/5;
    plusOneIndex = signal>threshold;
    minusOneIndex = signal<(-threshold);

    % Differences to detect 1 to -1 and visa versa or to 0. 
    diffSignal = [0, diff(signal)];
    DiffThresholdZero = max(abs(diffSignal))/2; 
    zeroIndex = signal<threshold & signal>-threshold & (diffSignal < DiffThresholdZero & diffSignal > -DiffThresholdZero);

    % Calculate Ternary bit sample hits
    lastInd = find(diffSignal<-DiffThresholdZero,1,'first');
    firstPeakID = lastInd - round((1/Ft)/2/objScope.sample_interval);
    ternaryTime=time(firstPeakID):(1/Ft):time(end);


    % decode to ternary bits
    ternaryBit = nan(1,length(signal));
    ternaryBit(plusOneIndex) = 0.1;
    ternaryBit(minusOneIndex) = -0.1;
    ternaryBit(zeroIndex) = 0;

    figure
    hold on
    plot(time(1:1000),signal(1:1000))
    plot(time(1:1000),signal(1:1000))
    plot(time(1:1000),ternaryBit(1:1000),'*');
    plot(time(1:1000),diffSignal(1:1000))
    plot(time(firstPeakID),signal(firstPeakID),'go')
    for i=1 : find(ternaryTime<= time(1000),1,'last')
        xline(ternaryTime(i));
    end
    yline(threshold,'--g','+Threshold')
    yline(-threshold,'--g','-Threshold')
end

function twoTernaryToDataSymbols(objScope)

end
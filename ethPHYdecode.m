classdef ethPHYdecode
    properties
        time;
        value;
    end
    methods
        function obj = ethPHYdecode(varargin)
           
            [objScope, chNr, threshold, cut_off_frequency, verbose] = obj.splitVarargin(varargin);
            
            %try
            if(verbose);disp(['Starting ETH decoding:']);tic;end;
            [X,Y] = obj.mlt(objScope,'channelnr',chNr,'threshold',threshold,'cut_off_frequency',cut_off_frequency,'verbose',verbose-(verbose>0));
            if(verbose);disp([' 1) MLT-decoding took ' mat2str(round(toc,3)) ' seconds.']);tic;end;
            [Y] = obj.descrambler(Y,verbose-(verbose>0));
            if(verbose);disp([' 2) Descrambling took ' mat2str(round(toc,3)) ' seconds.']);tic;end;
            [X,Y] = obj.decoder5B4B(X,Y,verbose-(verbose>0));
            if(verbose);disp([' 3) 5B to 4B decoding took ' mat2str(round(toc,3)) ' seconds.']);tic;end;
            %[X,Y] = convertPacketBytes(X,Y,verbose-(verbose>0));
            %if(verbose);disp([' 4) Nibble to Packet conversion took ' mat2str(round(toc,3)) ' seconds.']);tic;end;
            if(verbose);disp('ETH decoding finished succesfull.');end;
            obj.time = X;
            obj.value = Y;
            
            % catch ex % exception
            %     if strcmp(ex.identifier,'mlt:sampleRateToLow')
            %         warning off backtrace
            %         warning(ex.message);
            %         warning on backtrace
            %     else
            %         err.identifier = ex.identifier;
            %         err.message = ex.message;
            %         error(err);
            %     end
            % end
        end
    end
    methods (Static)
        function version
            disp('Version: 1.0');
            disp('Release Date: 1-03-2017');
        end
        function help
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            disp('%%                          ETHERNET DECODER                           %%');
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            disp('%                                                                       %');
            disp('% Decodes a 100Mbps ethernet signal from analog values to an eth object %');
            disp('%                                                                       %');
            disp('%  Internal Functions:                                                  %');
            disp('%    - <a href="matlab: help ethPHYdecode>mlt;">mlt</a>                                                              %');
            disp('%    - <a href="matlab: help ethPHYdecode>descrambler;">descrambler</a>                                                      %');
            disp('%    - <a href="matlab: help ethPHYdecode>decoder5B4B;">decoder5B4B</a>                                                      %');
            disp('%                                                                       %');
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        end
    end
    methods (Static, Hidden, Access=private)
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%                              FUNCTIONS                              %%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
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
        
      
        function [Vx,Vy] = mlt(varargin)
            [objScope, chNr, threshold, cut_off_frequency, verbose] = ethPHYdecode.splitVarargin(varargin);
            
            %% CHECK FOR A HIGH ENOUGH SAMPLE RATE
            if(objScope.sample_interval > 4e-9)
                errorStruct.message = 'Sample rate of at least 250 MHz is required';
                errorStruct.identifier = 'mlt:sampleRateToLow';
                error(errorStruct);
            end
            %% INTERPOLATE TILL 1GHZ SAMPLE RATE
            X = objScope.time;
            % For backward compatibility check if objscope contains value,
            % with new version the scope object contains multiple channels.
            if isprop(objScope,"value")
                Y = objScope.value(chNr);
            else
                Y = objScope.channels(chNr).value;
            end
            
            Scale = 1;
            mainClockSampleTime = 8e-9; % main Clock frequency of 100BASE-TX
            
            while objScope.sample_interval > 1e-9
                objScope.sample_length = objScope.sample_length * 2;
                objScope.sample_interval = objScope.sample_interval / 2;
                Scale = Scale * 2;
            end
            
            X = X(1)-objScope.sample_interval*(Scale-3)/2:objScope.sample_interval:X(end)+objScope.sample_interval*(Scale+1)/2;
            
            expfft = ceil(log2(length(Y)));
            fftResult = fft(Y,2^expfft);
            % Adding zero's (resulting in 1 GS/s)
            fftResult = [fftResult(1:floor((end+1)/2))*Scale, zeros(1,length(fftResult)*(Scale-1)), fftResult(ceil((end+1)/2):end)*Scale];
            Fs = 1/objScope.sample_interval; % Sampling frequency
            L = objScope.sample_length;      % Length of signal
            
            % Cutting 40-60Hz, 140-160Hz and Cut-off frequencies
            filterfreq = [40,60;140,160;cut_off_frequency,1e99];        
            for i = 1:size(filterfreq,1)
                locmin([1,2]) = round(filterfreq(i,1)*L/Fs)+1;
                locmax([1,2]) = round(filterfreq(i,2)*L/Fs)+1;
                if(locmin(1)<1);locmin(1)=1;end;
                if(locmax(1)<1);locmax(1)=1;end;
                if(locmin(2)<2);locmin(2)=2;end;
                if(locmax(2)<2);locmax(2)=2;end;
                if(locmin(1)>length(fftResult)/2);locmin(1)=round(length(fftResult)/2);end;
                if(locmax(1)>length(fftResult)/2);locmax(1)=round(length(fftResult)/2);end;
                if(locmin(2)>length(fftResult)/2);locmin(2)=round(length(fftResult)/2);end;
                if(locmax(2)>length(fftResult)/2);locmax(2)=round(length(fftResult)/2);end;
                L1 = locmin(1):locmax(1);
                L2 = length(fftResult)-(locmax(2)-2):length(fftResult)-(locmin(2)-2);
                fftResult(L1) = 0;
                fftResult(L2) = 0;
            end
            
            Yorig = Y;
            Y = real(ifft(fftResult,(2^expfft)*Scale));
            Y = Y(1:length(X));
            
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
                plt(X(plotMin:plotMax),Y(plotMin:plotMax),'downsample',1e5,'r');             
                xlabel('t(s)');
                ylabel('U(V)');
                axis([Xmin,Xmax,Ymin,Ymax]);
                hold off;
                
                figure;
                f = Fs*(0:ceil(L/2))/L;
                P1 = abs(fftResult(1:ceil(L/2)+1));
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
            clear F Yorig Scale filterfreq;
            
            %% PLOTTING GRAPH 1  
            if(verbose)
                figure;
                subplot(4,1,1);
                XT = objScope.time;
                
                if isprop(objScope,"value")
                    YT = objScope.value{chNr};
                else
                    YT = objScope.channels(chNr).value;
                end
                hold on;
                plot(XT(Xmin<XT & XT<Xmax),YT(Xmin<XT & XT<Xmax), 'color', [0.5 0.5 0.5]);
                plot(X(plotMin:plotMax),Y(plotMin:plotMax),'r');
                axis([Xmin,Xmax,Ymin,Ymax]);
                subtitle('Signals');
                legend('raw','upscaled&filtered');
                hold off;
            end
            objScope.time = X;
            
            if isprop(objScope,"value")
                objScope.value{chNr} = Y;
            else
                objScope.channels(chNr).value = Y;
            end
            
            %% DETECTING FALLING AND RISING EDGES
            X = X + objScope.sample_interval/2;
            Y = [abs(diff(Y)),0];
            Y = Y/max(Y);
            %% PLOTTING GRAPH 2-A
            if(verbose)
                subplot(4,1,2);
                hold on;
                plot(X(plotMin:plotMax),Y(plotMin:plotMax),'r');
                
            end
            %% CALCULATING ONE AND ZERO AREA'S
            
            %hl = floor((Freq/objScope.sample_interval)*0.9)+1;
            convWindow = [1 1 1 1 1 1 1 1];
            %convWindow = [0.2 0.5 0.75 1 1 0.75 0.5 0.2];
            
            
            convolutionOutput = conv(Y,convWindow,'same');
            sortConvOutput = sort(convolutionOutput);
            ConvOutMin = sortConvOutput(floor(length(sortConvOutput)*0.1));
            ConvOutmax = sortConvOutput(ceil(length(sortConvOutput)*0.9));
            convolutionOutput = (convolutionOutput-ConvOutMin)/(ConvOutmax-ConvOutMin);
            
            %Y = cY > 0.5;
            %cY = cY/max(cY);
            condition = true;
            loopcounter = 0;
            while(condition)
                bitLevel = convolutionOutput > threshold;
                meanBitLevel = mean(bitLevel);
                if((0.49 < meanBitLevel && meanBitLevel < 0.51) || loopcounter >= 0)
                    condition = false;
                else
                    loopcounter = loopcounter + 1;
                    threshold = threshold - 0.1*(0.5-meanBitLevel);
                end
            end
            if(verbose)
                disp(['Tuned threshold value: ' num2str(threshold)]);
                disp(['Resulting average: ' num2str(meanBitLevel)]);
            end
            clear meanY condition loopcounter             
            %% PLOTTING GRAPH 2-B
            if(verbose)
                plot(X(plotMin:plotMax),convolutionOutput(plotMin:plotMax),'color',[0.25 0.25 0.25]);
                plot(X(plotMin:plotMax),bitLevel(plotMin:plotMax),'b');
                axis([Xmin,Xmax,-0.25,1.25]);
                subtitle('Bit level calculation');
                legend('difference','conv out','Bit level');
                hold off;
            end
            %% DETERMINING RELIABLE BIT SAMPLES
            % detect changes
            dY = [diff(bitLevel),0];
            
            % shift time with half sample time to set the sample in the
            % middle
            shiftedTime = X + objScope.sample_interval/2;
            
            % get time on bit chanches
            timeRisingEdge = shiftedTime(dY(1:end-1)==+1);
            timeFallingEdge = shiftedTime(dY(1:end-1)==-1);
            % check length of time variables
            if length(timeFallingEdge) > length (timeRisingEdge); timeFallingEdge = timeFallingEdge(2:end); end
            if length(timeFallingEdge) < length (timeRisingEdge); timeRisingEdge = timeRisingEdge(2:end); end
            
            % get time of the centerpoint between rising and falling edge
            BitTiming = (timeRisingEdge + timeFallingEdge)/2;
            
%             figure
%             plot(bitLevel)
%             hold on
%             plot(dY)
%             legend('bitLevel','dY')
            
            dsX = (abs(timeRisingEdge-timeFallingEdge))/mainClockSampleTime;
            clear sXH sXL
            usX = BitTiming((1.995<dsX & dsX<2.005)); % EVEN
            usX = [usX - mainClockSampleTime/2 , usX + mainClockSampleTime/2];
            usX = [usX , BitTiming((0.995<dsX & dsX<1.005))]; % ODD
            BitTiming = sort(usX); % SORTING
            clear usX dsX;
            if(verbose);disp(['Number of reliable bit samples: ' num2str(length(shiftedTime))]);end
            
            
            %% DETERMINING UNRELIABLE BIT SAMPLES BY INTERPOLATION
            dsX = diff(BitTiming);
            j = 1;
            Vx = zeros(1,round((max(BitTiming)-min(BitTiming))/mainClockSampleTime+1));
            timeIndex = 1;
            for i = 1:1:length(BitTiming)-1
                timeIndex = round(dsX(i)/mainClockSampleTime);
                stepSize = dsX(i)/timeIndex;
                Vx(j:j+timeIndex) = BitTiming(i)+(0:stepSize:timeIndex*stepSize);
                j = j+timeIndex;
            end
            
            if length(Vx)>j+timeIndex; Vx(j+timeIndex+1:end) = []; end
            clear j timeIndex stepSize sX dsX;
            %% GETTING BIT VALUE FOR ALL DETERMINED SAMPLES
            Vy = zeros(1,length(Vx));
            j = find(X>Vx(1)-mainClockSampleTime/2,1);
            Interval = mainClockSampleTime/4;
            for i=1:length(Vx)
                k = 0;
                while X(j) < Vx(i) - Interval; j = j + 1; end
                while X(j) < Vx(i) + Interval
                    multiplier = 1-abs(Vx(i)-X(j)).^1.3/Interval;
                    Vy(i) = Vy(i) + multiplier*convolutionOutput(j);                 
                    k = k + multiplier;
                    j = j + 1;                    
                end
                Vy(i) = Vy(i)/k;
            end
            sVy = sort(Vy);
            sVy = sVy(floor(length(Vy)*0.1):ceil(length(Vy)*0.9));
            dsVy = sVy(2:end)-sVy(1:end-1);
            mdsVy = max(dsVy);
            fmdsVy = find((dsVy==mdsVy),1);
            VyThreshold = (sVy(fmdsVy)+sVy(fmdsVy+1))/2;
            clear X Y bitLevel k i j;            
            %% PLOTTING GRAPH 3
            if(verbose)
                subplot(4,1,3);
                hold on;
                Zeros = length(Vy(Vy<=0.4));
                Others = length(Vy(Vy>0.4 & Vy<0.6));
                Ones = length(Vy(Vy>=0.6));
                text(0.05,0.50,['Zeros: ' num2str(Zeros) '  (' num2str(100*Zeros/length(Vy)) '%)']);
                text(0.05,0.75,['Uncertain: ' num2str(Others)]);
                text(0.05,1.00,['Ones:  ' num2str(Ones) '  (' num2str(100*Ones/length(Vy)) '%)']);
                line([0 1],[VyThreshold VyThreshold],'Color',[0.5 0.5 0.5], 'LineWidth',1);
                plt((1:length(sVy))/length(sVy),sVy,'downsample',100000,'.b','MarkerSize',4);
                if(Others>0); plt((1:Others)/Others,sort(Vy(Vy>0.4 & Vy<0.6)),'.r','MarkerSize',4); end
                axis([0,1,-0.25,1.25]);
                hold off;
                subtitle('Sampled bit level');
%                 legend('difference','conv out','Bit level');
                clear Zeros Others Ones TotalBits;
            end
            %% DISCRETIZE BIT VALUES                     
            Vy = Vy > VyThreshold;
            %% PLOTTING GRAPH 4
            if(verbose)
                subplot(4,1,4);
                hold on;
                X = objScope.time;
                if isprop(objScope,"value")
                    Y = objScope.value{chNr};
                else
                    Y = objScope.channels(chNr).value;
                end
                plot(X(plotMin:plotMax),Y(plotMin:plotMax),'r');
                stem(Vx(Xmax>Vx & Vx>Xmin),Vy(Xmax>Vx & Vx>Xmin),'b');
                axis([Xmin,Xmax,Ymin,Ymax]);
                hold off;
                subtitle('Bit level vs signal');
                legend('signal','bit level');
                clear X Y ;
            end
        end
        %%
        function [plainText] = descrambler(cipherText,verbose)
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %%                             DESCRAMBLER                             %%
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %                                                                       %
            % % Descrambles a pseudo-random scrabled bit matrix                       %
            % % using a 11 bit shift register                                         %
            % %                                                                       %
            % %  Parent Function: <a href="matlab: help ethPHYdecode;">ethPHYdecode</a>                                         %
            % %                                                                       %
            % %                                                                       %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�%%%%
            if(~exist('verbose','var'));verbose=-1;warn('All underlying functions are executed in verbose mode');end;
            finished=0;
            i=0;
            sync=0;
            while finished==0
                % FIND SYNC PATTERN
                j = 11+i;
                keyStream = cipherText(j-10:j);
                descrambledKey = [0 0 0 0 0 0 0 0 0 0 0];
                while sync==0
                    if j == length(cipherText)
                        sync = 0;
                        finished = 1;
                        break;
                    end
                    j = j+1;
                    keyStream(2:12) = keyStream(1:11);
                    keyStream(1) = cipherText(j);
                    descrambledKey(2:11) = descrambledKey(1:10);
                    % calculate result of 1+ x^9+x^11
                    descrambledKey(1) = xor(keyStream(1), xor(keyStream(10), keyStream(12)));
                    
                    % PATTERN MATCH
                    if descrambledKey == [1 1 1 1 1 1 1 1 1 1 1]
                        keyStream = ~keyStream(1:11);
                        syncLoc = j;
                        syncK = keyStream;
                        sync=1;
                    end
                end
                i = i+1;
                
                % DESCRAMBLE
                plainText = zeros(1,length(cipherText));
                for j = syncLoc:length(cipherText)
                    keyStream(2:12) = keyStream(1:11);
                    keyStream(1) = xor(keyStream(12), keyStream(10));
                    plainText(j) = xor(cipherText(j), keyStream(1));
                    
                    if j > 100
                        finished = 1;
                        break;
                    elseif plainText(j) == 1
                        %IDLE = OK
                    elseif plainText(j-6:j) == [1 1 1 1 1 1 0]
                        %J ?
                    elseif plainText(j-7:j) == [1 1 1 1 1 1 0 0]
                        %J ?
                    elseif plainText(j-8:j) == [1 1 1 1 1 1 0 0 0]
                        %K?
                    elseif plainText(j-10:j) == [1 1 1 1 1 1 0 0 0 1 0]
                        %K?
                    elseif plainText(j-11:j) == [1 1 1 1 1 1 0 0 0 1 0 0]
                        %K?
                    elseif plainText(j-12:j) == [1 1 1 1 1 1 0 0 0 1 0 0 0]
                        %K?
                    elseif plainText(j-13:j) == [1 1 1 1 1 1 0 0 0 1 0 0 0 1]
                        %K = OK
                        finished = 1;
                    else
                        sync = 0;
                        break;
                    end
                end
            end
            % CALCULATE BACK TO 1
            keyStream = syncK;
            for a=1:syncLoc
                keyStream(12) = xor(keyStream(10), keyStream(1));
                keyStream(1:11) = keyStream(2:12);
            end
            
            plainText = zeros(1,length(cipherText));
            % DESCRAMBLE ALL
            for i = 1:length(cipherText)
                for j = 12:-1:2
                    keyStream(j)=keyStream(j-1);
                end
                keyStream(1) = mod(keyStream(12)+keyStream(10),2);
                plainText(i) = mod(cipherText(i)+keyStream(1),2);
            end
        end
        %%
        function [Vx,Vy] = decoder5B4B(X,Y,verbose)
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %%                          5B to 4B DECODER                           %%
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %                                                                       %
            % % Decodes 5B groups to 4B nibbles                                       %
            % %                                                                       %
            % %  Parent Function: <a href="matlab: help ethPHYdecode;">ethPHYdecode</a>                                         %
            % %                                                                       %
            % %                                                                       %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�%%%%
            warning off backtrace;
            if(~exist('verbose','var'));verbose=-1;warning('All underlying functions are executed in verbose mode');end;
            loc = find(~Y);
            loc = loc(loc>12 & loc<(length(Y)-7));
            align = -1;
            for i = 1:length(loc)
                if Y(loc(i)-12:loc(i)+7) == [1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 1 0 0 0 1]
                    align = mod(loc(i)-3,5)+1;
                    break;
                end
            end
            if align == -1
                warning('No packets found');
                align = 1;
            end
            Vx = zeros(1,floor((length(Y)-align)/5));
            Vy = zeros(1,floor((length(Y)-align)/5));
            %%% Decoding 5B to 4B
            for i = 1:floor((length(Y)-align)/5)
                Vx(i) = X(align+(i-1)*5);
                Vy(i) = ethPHYdecode.decoder5B4B_LUT(Y(align+(i-1)*5:align+(i-1)*5+4));
            end
            %%% Detecting Packets
            if(verbose)
                disp(['H= ' num2str(length(find(Vy==17)))]);
                disp(['I= ' num2str(length(find(Vy==18)))]);
                disp(['J= ' num2str(length(find(Vy==-1)))]);
                disp(['K= ' num2str(length(find(Vy==-2)))]);
                disp(['Q= ' num2str(length(find(Vy==26)))]);
                disp(['R= ' num2str(length(find(Vy==-4)))]);
                disp(['S= ' num2str(length(find(Vy==28)))]);
                disp(['T= ' num2str(length(find(Vy==-3)))]);
                disp(['Nan= ' num2str(length(find(Vy==32)))]);
            end
        end
        %%
        function [D] = decoder5B4B_LUT(E)
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %%                 5B to 4B DECODER Look Up TABLE                      %%
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            switch(E(1)*16+E(2)*8+E(3)*4+E(4)*2+E(5))
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 31; D = 18; % Idle (I)
                case 4;  D = 17; % Transmit Error (H)
                case 24; D = -1; % SSD (Part 1) (J)
                case 17; D = -2; % SSD (Part 2) (K)
                case 0;  D = 26; % (Q)
                case 7;  D = -4; % ESD (Part 2) (R)
                case 25; D = 28; % (S)
                case 13; D = -3; % ESD (Part 1) (T)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 30; D = 0;  % Data 0
                case 9;  D = 1;  % Data 1
                case 20; D = 2;  % Data 2
                case 21; D = 3;  % Data 3
                case 10; D = 4;  % Data 4
                case 11; D = 5;  % Data 5
                case 14; D = 6;  % Data 6
                case 15; D = 7;  % Data 7
                case 18; D = 8;  % Data 8
                case 19; D = 9;  % Data 9
                case 22; D = 10; % Data A
                case 23; D = 11; % Data B
                case 26; D = 12; % Data C
                case 27; D = 13; % Data D
                case 28; D = 14; % Data E
                case 29; D = 15; % Data F
                otherwise; D = 32; % NOT IN THE LIST (NaN)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
        
        %%
        function [Vx,Vy] = decoder4B3B(X,Y,verbose)
            
            
        end
        %%
        function [Vx,Vy] = convertPacketBytes(X,Y,verbose)
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %%                Nibble to Bytes per Packet converter                 %%
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %                                                                       %
            % % Converts the Nibbles into Bytes                                       %
            % %                                                                       %
            % %  Parent Function: <a href="matlab: help ethPHYdecode;">ethPHYdecode</a>                                         %
            % %                                                                       %
            % %                                                                       %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�%%%%
            warning off backtrace;
            if(~exist('verbose','var'));verbose=-1;warning('All underlying functions are executed in verbose mode');end;
            X = X(Y<16);
            Y = Y(Y<16);
            J = find(Y==-1);
            K = find(Y==-2);
            R = find(Y==-3);
            S = find(Y==-4);
            if range(mod([J,R,K+1,S+1],2)) == 0
                Index = mod(J(1)-1,2)+1;
                Vx = X(Index:2:end);
                Vy = Y(Index:2:end)*16+Y(Index+1:2:end);
            else
                error('Packets need to be determined seperately');
            end
        end
    end
end

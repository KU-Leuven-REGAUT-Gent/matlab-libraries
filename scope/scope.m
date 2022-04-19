%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%                           SCOPE CLASS                              %%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %                                                                      %
%  %    Author: Frederic Depuydt                                          %
%  %    Company: KU Leuven                                                %
%  %    Contact: frederic.depuydt@kuleuven.be; f.depuydt@outlook.com      %
%  %    Version: 1.3                                                      %
%  %                                                                      %
%  %    An Scope class to analyse scope signals                           %
%  %    Readable files: isf, csv                                          %
%  %                                                                      %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %        FUNCTIONS (static)                 *Object creation*          %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %                                                                      %
%  %    usage: objScope = scope.function(var1, var2, ...)                 %
%  %                                                                      %
%  %    isfread(                Reading scope signals from an ISF file    %
%  %        file,                   Filename + extension as String        %
%  %        verbose)                Integer to enable verbose mode        %
%  %                                                                      %
%  %    wfmread(                Reading scope signals from an WFM file    %
%  %        file,                   Filename + extension as String or directory      %
%  %        verbose)                Integer to enable verbose mode        %
%  %                                       option to specify measurement and extract all corresponding channels                               %
%  %    csvread(                DEPRECATED! Reading from a CSV file       %
%  %        file,                   Filename + extension as String        %
%  %        channels,               Array of strings refering to channels %
%  %        verbose,                Integer to enable verbose mode        %
%  %        retime)                 Calculate more accurate timestamps    %
%  %                                                                      %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %        FUNCTIONS (non-static)                                        %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %                                                                      %
%  %    usage: result = objScope.function(var1, var2, ...)                %
%  %    note: most functions do not alter the original scope object,      %
%  %          but return a new object with the function results           %
%  %                                                                      %
%  %    values(                 Returning the values of a channel         %
%  %        channels)               Array of strings refering to channels %
%  %            returns: matrix of the requested values                   %
%  %                                                                      %
%  %    fft(                 Returning the fft of a channel  (with zero padding)      %
%  %        scale,               'db' converts amplitude to decibels %
%  %        window,              hanning, hamming, blackmanharris, gaussian, rectangle, none is possible  %
%  %            returns:  [frequency axis, spectrum, window size]                  %
%  %                                                                      %
%  %    plotChannels(                 plot the scope signals
%  %        'channels', channelArray              array of channels to plot
%  %        'titles', titleArray                     array of subtitles Standard CHx - title e.g. ["test";"test2"]
%  %        'limit', xlimitArray                   array with start point and end point of x axis
%  %        "save"  , saveName                    when the string save is detected the next attribute  (if it's a string) will be
%  %                                        the filename of the stored figures with _channels after it.
%  %                                        Otherwise the filename of the scope will be used.
%  %            e.g     s.plotChannels('channesls,[2:3], 'save','test','titles'["title for channel 2"; "title for channel 3"])
%  %                    %
%  %    plotMath(                 plot the scope signals
%  %        'channels', channelArray              array of channels to plot
%  %        'titles', titleArray                    array of subtitles Standard CHx - title  e.g. ["test";"test2"]
%  %        'limit', xlimitArray                    array with start point and end point of x axis
%  %        "save"  , saveName                  when the string save is detected the next attribute  (if it's a string) will be
%  %                                        the filename of the stored figures with _math after it.
%  %                                        Otherwise the filename of the scope will be used.
%  %            e.g     s.plotMath('channesls,[2:3], 'save','test','titles'["title for channel 2"; "title for channel 3"])
%  %                   png and fig are saved
%  %                                                                      %
%  %    split(                  Splitting 1 scope object into 2 (a and b) %
%  %        channels_a,             Array of strings refering to channels %
%  %        channels_b)             Array of strings refering to channels %
%  %            returns: [objScope1, objScope2]                           %
%  %                                                                      %
%  %    downsample(             Lowering the amount of sample rate        %
%  %        samples,                The new number of samples             %
%  %        verbose)                Integer to enable verbose mode        %
%  %            returns: downsampled scope object                         %
%  %                                                                      %
%  %    remove(                 Removing channels from a Scope object     %
%  %        channels,               Array of strings refering to channels %
%  %        verbose)                Integer to enable verbose mode        %
%  %            returns: scope object without the removed channels        %
%  %                                                                      %
%  %    noisefilter(            Filters noise from requested channels     %
%  %        values,                 The input values you want to filter   %
%  %        threshold,              Threshold value to be filtered        %
%  %        verbose)                Integer to enable verbose mode        %
%  %            returns: filtered output values                           %
%  %                                                                      %
%  %    bandstop(               Filtering by frequencybands               %
%  %        values,                 The input values you want to filter   %
%  %        freq,                   Array of frequentybands to be filtered%
%  %        verbose)                Integer to enable verbose mode        %
%  %            returns: filtered output values                           %
%  %                                                                      %
%  %    scale(                  Scaling values of the requested channels  %
%  %        channels,               Array of strings refering to channels %
%  %        target,                 Target value to which to scale to     %
%  %        verbose)                Integer to enable verbose mode        %
%  %            returns: nothing, results directly applied on scope object%
%  %                                                                      %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %                                                                      %
%  %    VERBOSE MODE: (default=-1)                                        %
%  %        all static functions check for a verbose variable             %
%  %        to enable or disable display output to the console            %
%  %                                                                      %
%  %    verbose ==  0;  % Display output disabled                         %
%  %    verbose == -1;  % Display output enabled for all nested functions %
%  %    verbose ==  x;  % Display output enabled for x nested functions   %
%  %                                                                      %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef scope < dynamicprops & matlab.mixin.Copyable
    properties
        model, ...
            fileName,...
            firmware_version, ...
            waveform_type, ...
            point_format, ...
            horizontal_units, ...
            horizontal_scale, ...
            horizontal_delay, ...
            sample_interval,...
            filter_frequency, ...
            record_length, ...
            sample_length, ...
            gating, ...
            gating_min, ...
            gating_max, ...
            probe_attenuation, ...
            vertical_units, ...
            vertical_offset, ...
            vertical_scale, ...
            vertical_position, ...
            time, ...
            channels, ...
            
    end
    methods
        function obj = scope(model)
            obj.model = model;
        end
        function result = values(obj,str)
            for i=1:numel(obj.channels)
                if(strcmp(str,obj.channels{i}.name))
                    result = obj.channels{i}.value;
                end
            end
        end
        % ----------------------- FFT function ----------------------------
        function [freq_axis, fft_result,N] = advancedFFT(obj,scale,window,gatePosition,gateDuration,levelOffset)
            
            recordStart = gatePosition-gateDuration/2;
            recordEnd = gatePosition+gateDuration/2;
            periodExtracted = zeros(1,length(obj.time));
            t1 = find(obj.time>=recordStart,1,'first');
            t2= find(obj.time>=recordEnd,1,'first');
            periodExtracted(t1:t2) = 1;
            figure
            hold on
            plot(objScope(1).time,obj.channels{1}.value)
            plot(objScop.time,periodExtracted);
            hold off
            
            sizeData = length(objScope(1).value{1}(t1:t2-1));
            N=sizeData;%2^(nextpow2( sizeData));
            Fs=1/objScope(1).sample_interval;
            
            if sizeData< N
                data =zeros(1,N);
                data(1:sizeData)= obj.channels{1}.value(t1:t2-1);
            else
                data= obj.channels{1}.value(t1:t2-1);
                data=data(1:N) ;
            end
            
            switch window
                case 'hanning'
                    %Characteristics:
                    % -     Better frequency,
                    % -     Poorer magnitude accuracy than Rectangular.
                    % -     Hanning has slightly poorer frequencyresolution than Hamming.
                    % Best for:
                    % -     Sine, periodic, and narrow-band random noise.
                    % -     Transients or bursts where the signal levels before and after the event are significantly different
                    fftWindow = hanning(N);
                case 'hamming'
                    %Characteristics:
                    % -     Better frequency, ,
                    % -     poorer magnitude accuracy than Rectangular
                    % -     THamming has slightly better frequencyresolution than Hanning
                    % Best for:
                    % -     Sine, periodic, and narrow-band random noise
                    % -     Transients or bursts where the signal levels before and after the event are significantly different
                    fftWindow = hamming(N);
                case 'rectangle'
                    %Characteristics:
                    % -     Best frequency,
                    % -     worst magnitude resolution
                    % -     This isessentially the same as no window
                    % Best for:
                    % -     Transients or bursts where the signal levels before and after the event are nearly equal
                    % -     Equal-amplitude sine waves with frequencies that are very close
                    fftWindow = rectwin(N);
                case 'blackmanharris'
                    %Characteristics:
                    % -     worst frequency resolution,
                    % -     Best magnitude
                    % -     This isessentially the same as no window
                    % Best for:
                    % -     Predominantly single frequency signals to look forhigher order harmonics
                    fftWindow = blackmanharris(N);
                case 'gaussian'
                    %Characteristics:
                    % -     worst frequency resolution,
                    % -     Best magnitude
                    % -     This isessentially the same as no window
                    % Best for:
                    % -     Predominantly single frequency signals to look forhigher order harmonics
                    fftWindow = gausswin(N);
                otherwise
                    
            end
            if exist('fftWindow','var')
                data = data'.*fftWindow;
            else
                data = data';
            end
            Fs=1/obj(1).sample_interval;
            spectr_res= Fs/N;
            max_freq_possible = Fs/2;
            freq_axis = (Fs*(0:(N/2))/N)';
            
            fft_res = fft(data,N);
            fft_abs= (abs(fft_res/N));
            fft_result=fft_abs(1:N/2+1);
            fft_result(2:end-1) = 2*fft_result(2:end-1);
            if strcmp(scale, 'db')
                fft_result = mag2db(fft_result)- mag2db(levelOffset);
            end
            
        end
        
        function fft(obj)
            figure;
            fontSize = 20;
            set(gca,'fontsize',fontSize+2) % set fontsize of the plot to 20
            set(gcf,'units','normalized','outerposition',[0 0 1 1]) % full screen
            set(0, 'DefaultAxesFontSize', fontSize);
            subplotArray=[];
            chSize = numel(obj.channels);
           
            [freq_axis, fft_result] = obj.channels(1).fft(obj.sample_interval,"","hamming");
            subplotArray(1) =subplot(chSize,1,1);
            plot(freq_axis, fft_result,'LineWidth',2)
             title('FFT of scope channels');
            subtitle(['Channel ' obj.channels(1).name])
            
            ylabel('amplitude');
            for i=2:chSize
                [freq_axis, fft_result] = obj.channels(i).fft(obj.sample_interval,"db","hamming");
                subplotArray(i) =subplot(chSize,1,i);                
                plot(freq_axis, fft_result,'LineWidth',2)
                subtitle(['Channel ' obj.channels(i).name])
                ylabel('amplitude');
            end
            xlabel('Frequency [Hz]');
            
            linkaxes(subplotArray,'x');
        end
        
        % ----------------------- plot function ----------------------------
        
        function plotChannels(obj,varargin)
            % input arguments
            % 1 - single value or array of channels that needs to be
            % plotted
            % 2 - title array ( default CHx)
            % 3 - y label array (default Voltage [V]
            %------------------------------------------------------
            
            % Declaration and initialisation of the titles, ylabels and
            % channels variables
            [ch, xLimits, savePlot, saveName, titles] = scope.splitVarargin(varargin);
            if isempty(ch)
                ch = str2double(extractAfter({obj.channels.name},'Ch')');
            end
            saveName = strcat(saveName,'_channels');
            chSize = numel(ch);
            chScope= figure('name',['scope plot measurement ' obj.fileName ]);
            fontSize = 20;
            set(gca,'fontsize',fontSize+2) % set fontsize of the plot to 20
            set(gcf,'units','normalized','outerposition',[0 0 1 1]) % full screen
            set(0, 'DefaultAxesFontSize', fontSize);
            subplotArray=[];
            
            s=1;
            for i=1:numel(ch)
                for  j =1:numel(obj.channels)
                    if contains(obj.channels(j).name,string(ch(i)))
                        % declaration ylabel
                        switch obj.channels(j).vertical_unit
                            case 'V'
                                yText = ["Voltage [" + obj.channels(j).vertical_unit  + "]"] ;
                            case 'A'
                                yText = ["Current [" + obj.channels(j).vertical_unit  + "]"] ;
                        end
                        
                        subplotArray(s) =subplot(chSize,1,s);
                        plot(obj.time,obj.channels(j).value,'LineWidth',2,'Color',pltColor(ch(i)))
                        
                        ylabel(yText)
                        
                        if i<=numel(titles) && titles(i) ~= ""
                            title(obj.channels(j).name + " - " +  titles(i));
                        else
                            title(obj.channels(j).name);
                        end
                        
                        if ~isempty(xLimits)
                            xlim(xLimits);
                        else
                            xlim([obj.time(1), obj.time(end)]);
                        end
                        
                        s=s+1;
                        break;
                    end
                end
            end
            xlabel(["time ["+ obj.horizontal_units + "]"]);
            linkaxes(subplotArray,'x');
            %------- save plot ---------
            if savePlot
                D = pwd;
                if ~exist([D '\matlab'], 'dir')
                    mkdir([D '\matlab'])
                end
                
                print(chScope,'-dpng',fullfile(D,'matlab', saveName),'-r400');
                saveas(chScope,fullfile(D,'matlab', strcat(saveName, ".fig")));
            end
        end
        
        function plotMath(obj,varargin)
            % input arguments
            % 1 - single value or array of channels that needs to be
            % plotted
            % 2 - title array ( default CHx)
            % 3 - y label array (default Voltage [V]
            %------------------------------------------------------
            
            % Declaration and initialisation of the titles, ylabels and
            % channels variables
            [ch, xLimits, savePlot, saveName, titles] = scope.splitVarargin(varargin);
            saveName = strcat(saveName,'math');
            chSize = numel(ch);
            
            mathScope= figure('name',['scope math plot measurement ' obj.fileName ]);
            fontSize = 20;
            set(gca,'fontsize',fontSize+2) % set fontsize of the plot to 20
            set(gcf,'units','normalized','outerposition',[0 0 1 1]) % full screen
            set(0, 'DefaultAxesFontSize', fontSize);
            subplotArray=[];
            
            s=1;
            for i=1:numel(ch)
                for  j =1:numel(obj.math.channel)
                    if contains(obj.math(j).channel(j).name,string(ch(i)))
                        % declaration ylabel
                        switch obj.math(j).channel.vertical_unit
                            case 'dB'
                                yText =[ "Amplitude [" +  obj.math(j).channel.vertical_unit + "]"] ;
                                xText = ["frequency ["+ obj.math(j).horizontal_units + "]"];
                            otherwise
                                yText =obj.math(j).channel.vertical_unit ;
                                xText = obj.math(j).horizontal_units;
                        end
                        
                        subplotArray(s) =subplot(chSize,1,s);
                        plot(obj.math(j).time,obj.math(j).channel.value,'LineWidth',2)
                        
                        ylabel(yText)
                        
                        if i<=numel(titles) && titles(i) ~= ""
                            title(obj.math(j).channel.name + " - " +  titles(i,:));
                        else
                            title(obj.math(j).channel.name);
                        end
                        
                        if numel(ch) >1
                            
                            for  c=1:numel(ch)
                                Hunits(c,:) =  obj.math(c).horizontal_units;
                            end
                            % show only the xlabel on the lowest subplot
                            if sum(contains(cellstr(Hunits),Hunits(1,:)))==numel(ch)
                                if j == numel(ch) % show xlabel on last plot
                                    xlabel(xText)
                                end
                            else % show ylabel on each subplot when they differs
                                xlabel(xText)
                            end
                        else % only one plot
                            xlabel(xText)
                        end
                        
                        if ~isempty(xLimits)
                            xlim(xLimits);
                        else
                            xlim([obj.math(j).time(1), obj.math(j).time(end)]);
                        end
                        s=s+1;
                        break;
                    end
                end
            end
            linkaxes(subplotArray,'x');
            %------- save plot ---------
            if savePlot
                D = pwd;
                if ~exist([D '\matlab'], 'dir')
                    mkdir([D '\matlab'])
                end
                
                print(mathScope,'-dpng',fullfile(D,'matlab', saveName),'-r400');
                saveas(mathScope,fullfile(D,'matlab', strcat(saveName, ".fig")));
            end
        end
        
        %          function obj = copy(obj,parent)
        %              obj = parent;
        %          end
        
        
        % ----------------------- PROFINET decode function ----------------------------
        function obj = decodePN(obj,varargin)
            % Declaration and initialisation of channels
            switch nargin
                case 1
                    ch = 1:4;
                case 2
                    ch = varargin{1};
            end
            verbose =0;
            % decode the signals in the same order as the ch array.
            for i =1:numel(ch)
                for j=1: length(obj.channels)
                    if contains(obj.channels(j).name,string(ch(i)))
                        %              obj.pn.(obj.channels(j).name)  = eth.scoperead(obj,j,verbose);
                        obj.channels(j) = obj.channels(j).decodeChannelPN(obj,j,verbose);
                    end
                end
            end
        end
        
        function [obj1,obj2] = split(obj,channels1,channels2)
            obj1 = obj;
            obj2 = obj;
            obj1.channels.name = {};
            obj2.channels.name = {};
            obj1.channels.value = {};
            obj2.channels.value = {};
            k1 = 1;
            k2 = 1;
            for i=1:length(obj.channels)
                for j=1:length(channels1)
                    if(strcmp(channels1(j),obj.channels{i}))
                        obj1.channels{k1}.name=obj.channels{i}.name;
                        obj1.channels{k1}.value=obj.channels{i}.value;
                        k1 = k1+1;
                    end
                end
                for j=1:length(channels2)
                    if(strcmp(channels2(j),obj.channels{i}))
                        obj2.channels{k2}.name=obj.channels{i}.name;
                        obj2.channels{k2}.value=obj.channels{i}.value;
                        k2 = k2+1;
                    end
                end
            end
        end
        function obj = downsample(obj,samples,verbose)
            if(~exist('verbose','var'));verbose=-1;warn('All underlying functions are executed in verbose mode');end;
            if(~exist('samples','var'));return;end;
            n = obj.sample_length/samples;
            if(n>1 && mod(n,1)==0)
                newobj = obj;
                newobj.sample_interval = n*obj.sample_interval;
                newobj.sample_length = obj.sample_length/n;
                obj.time = mean(reshape(obj.time,n,[]));
                for i=1:length(obj.channels)
                    obj.channels{i}.value =  mean(reshape(obj.channels{i}.value,n,[]));
                end
            else
                error('Impossible downsample');
            end
        end
        
        function obj = remove(obj,channels,verbose)
            if(~exist('verbose','var'));verbose=-1;warn('All underlying functions are executed in verbose mode');end;
            if(~exist('channels','var'));return;end;
            for j=1:length(channels)
                for i=length(obj.channels):-1:1
                    if(strcmp(channels(j),obj.channels(i).name))
                        obj.channels(i).name=[];
                        obj.channels(i).value=[];
                    end
                end
            end
        end
        
        function Y = noisefilter(obj,X,threshold,verbose)
            if(~exist('verbose','var'));verbose=-1;warn('All underlying functions are executed in verbose mode');end;
            if(~exist('threshold','var'));threshold=0.0001;end;
            % Initial Values
            Fs = 1/obj.sample_interval;  % Sampling frequency
            T = 1/Fs;                    % Sampling period
            L = obj.sample_length;       % Length of signal
            F = fft(X);
            P1 = abs(F/L);
            F(P1<threshold)=0;
            Y = ifft(F);
            if(verbose)
                figure;
                P2 = abs(fft(Y)/L);
                ax = plt.getaxis(obj.time,[X,Y]);
                subplot(2,1,1);
                plt(obj.time,X,'downsample',1e6,color.ch1);
                title('Original signal');
                xlabel('t(s)');
                ylabel('X(t)');
                axis(ax);
                subplot(2,1,2);
                plt(obj.time,Y,'downsample',1e6,color.ch1);
                title('Filtered signal');
                xlabel('t(s)');
                ylabel('X(t)');
                axis(ax);
                
                figure;
                f = Fs*(0:ceil(L/2))/L;
                subplot(2,1,1);
                P1 = P1(1:ceil(L/2)+1);
                ax = [1 max(f) min(P1(P1>10^(-10))) 1];
                loglog(f,P1,color.ch2);
                title('Single-Sided Amplitude Spectrum of P1(t)');
                xlabel('f (Hz)');
                ylabel('|P1(f)|');
                axis(ax);
                subplot(2,1,2);
                P2 = P2(1:ceil(L/2)+1);
                loglog(f,P2,color.ch2);
                title('Single-Sided Amplitude Spectrum of P2(t)');
                xlabel('f (Hz)');
                ylabel('|P2(f)|');
                axis(ax);
            end
        end
        
        function Y = bandstop(obj,X,freq,verbose)
            if(~exist('verbose','var'));verbose=-1;warn('All underlying functions are executed in verbose mode');end;
            if(verbose);tic;end;
            if(~exist('freq','var')); Y=X; warn('No frequency range given, returned without filtering');return;end;
            if(size(freq,2) ~= 2); Y=X; warn('Frequency range invalid, returned without filtering');return;end;
            % Initial Values
            Fs = 1/obj.sample_interval;  % Sampling frequency
            L = obj.sample_length;       % Length of signal
            F = fft(X);
            P1 = F;
            for i = 1:size(freq,1)
                locmin([1,2]) = round(freq(i,1)*L/Fs)+1;
                locmax([1,2]) = round(freq(i,2)*L/Fs)+1;
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
            Y = real(ifft(F));
            if(verbose)
                toc
                P1 = abs(P1/L);
                P2 = abs(fft(Y)/L);
                figure;
                %ax = plt(obj.time,[X,Y]);
                subplot(2,1,1);
                plt(obj.time,X,'downsample',1e5,color.ch1);
                title('Original signal');
                xlabel('t(s)');
                ylabel('X(t)');
                %axis(ax);
                subplot(2,1,2);
                plt(obj.time,Y,'downsample',1e5,color.ch1);
                title('Filtered signal');
                xlabel('t(s)');
                ylabel('X(t)');
                %axis(ax);
                
                figure;
                f = Fs*(0:ceil(L/2))/L;
                subplot(2,1,1);
                P1 = P1(1:ceil(L/2)+1);
                ax = [1 max(f) min(P1(P1>10^(-10))) 1];
                loglog(f,P1,color.ch2);
                title('Single-Sided Amplitude Spectrum of P1(t)');
                xlabel('f (Hz)');
                ylabel('|P1(f)|');
                axis(ax);
                subplot(2,1,2);
                P2 = P2(1:ceil(L/2)+1);
                loglog(f,P1,color.lightgrey);
                axis(ax);
                hold on;
                loglog(f,P2,color.ch2);
                title('Single-Sided Amplitude Spectrum of P2(t)');
                xlabel('f (Hz)');
                ylabel('|P2(f)|');
                axis(ax);
                hold off;
            end
        end
        
        function Y = getValues(obj,channel)
            if(~exist('channel','var'));error('No channel selected');end;
            if(isnumeric(channel))
                Y = obj.channels(channel).value;
            else
                Y = obj.channels(channel).value;
            end
            
        end
        
        function scale(obj,str,target)
            if(~exist('pass','var'));target=0;end;
            for i=1:length(obj.channels)
                if(strcmp(str,obj.channels{i}.name))
                    YS = sort(obj.channels{i}.name,1, 'ascend');
                    fault = mean(YS(1:ceil(0.03*size(YS,1))))-target;
                    obj.channels{i}.value = obj.channels{i}.value - fault;
                end
            end
        end
    end
    methods (Access = protected)
        function thiscopy = copyElement(this)
            thiscopy = copyElement@matlab.mixin.Copyable(this); %shallow copy of all elements
            thiscopy.channels = copy(this.channels); %Deep copy of channels
        end
    end
    methods (Static)
        function [ch, xLimits, savePlot, saveName, titles] = splitVarargin(varargin)          
            varargin = varargin{1};
            if(numel(varargin) >= 1)
                while ~isempty(varargin)
                    if isa(varargin{1}, 'double') && ~exist('ch')
                        ch =varargin{1};
                        varargin(1) = [];
                    elseif isempty(varargin{1})
                        varargin(1) = [];
                    elseif(ischar(varargin{1}))
                        switch lower(varargin{1})
                            case 'channels'
                                ch = varargin{2};
                                varargin(1:2) = [];
                            case 'limit'
                                xLimits = varargin{2};
                                varargin(1:2) = [];
                            case 'save'
                                saveName = varargin{2};
                                savePlot = true;
                                varargin(1:2) = [];
                            case 'titles'
                                titles = varargin{2};
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
            
            if(~exist('ch','var'))
                ch = [];
            end
            if(~exist('saveName','var'))
                saveName = [];
                savePlot = false;
            end
            if(~exist('titles','var'))
                titles = [];
            end
            if(~exist('xLimits','var'));xLimits= [];end
        end
        
        function obj = isfread(file, verbose,varargin)
            if(~exist('verbose','var'));verbose=-1;warn('All underlying functions are executed in verbose mode');end;
            if ~exist('file', 'var')
                error('No file name, directory or pattern was specified.');
            end
            sizeRowVarargin = cellfun(@(c) size( c,1), varargin, 'UniformOutput', false);
            appIndex= find(cellfun(@(x)~isempty(strfind("app",x)), varargin( [sizeRowVarargin{:}]==1)));
            % Check whether file is a folder.
            if (isa(file, 'char' ) && exist(file, 'dir') )
                folder = file;
                % Get a list of all files that have the extension '.isf' or '.ISF'.
                files = [ dir(fullfile(folder, '*.isf')) ];
            elseif  isempty(appIndex) || appIndex==0
                % The pattern is not a folder, so must be a file name or a pattern with
                % wildcards (such as 'Traces/TEK0*.ISF').
                [folder, ~, ~] = fileparts(file);
                % Get a list of all files and folders which match the pattern...
                filesAndFolders = dir(file);
                % ...then exclude the folders, to get just a list of files.
                files = filesAndFolders(~[filesAndFolders.isdir]);
            end
            if (isempty(appIndex) || appIndex==0)
                fileNames = {files.name};
                datetimes = datestr([files.datenum]);
            end
            obj = scope('Unknown (ISF)');
            obj.firmware_version = 'Unknown (ISF)';
            if appIndex>0
                folder = varargin{appIndex +1};
                
                
                chNames = file(~cellfun('isempty',(regexpi(file,'CH'))));
                mthNames = file(~cellfun('isempty',(regexpi(file,'MTH'))));
                
                obj.fileName = extractBefore(file{1},'CH');
            elseif numel(fileNames)==0
                error('The pattern did not match any file or files: %s', file);
            elseif numel(fileNames) >1
                if contains(input('Extract one measurement: ','s'),["y","Y","yes","j","ja"])
                    fileInDir(:,2) = fileNames';
                    fileInDir(:,1)= num2cell(1:numel(fileNames))';
                    table(fileInDir)
                    fileNr =input('Measurement ID: ');
                    fName = extractBefore(fileInDir(fileNr,2),'CH');
                    fileIDs=~cellfun('isempty',regexp(fileNames,fName));
                    fileNames= fileNames(fileIDs);
                    chNames = fileNames(~cellfun('isempty',(regexp(fileNames,'tek[0-9]*CH'))));
                    mthNames = fileNames(~cellfun('isempty',(regexp(fileNames,'tek[0-9]*MTH'))));
                    obj.fileName = fName{1};
                else
                    warning('Returning empty scope object')
                    return;
                end
            elseif numel(fileNames) ==1
                fName = extractBefore(fileNames,'CH');
                if isempty(fName)
                   fName = extractBefore(fileNames,'MTH'); 
                end
                chNames = fileNames(~cellfun('isempty',(regexp(fileNames,'tek[0-9]*CH'))))
                mthNames = fileNames(~cellfun('isempty',(regexp(fileNames,'tek[0-9]*MTH'))))
                obj.fileName = fName{1};
                
            end
            
            % read setup file 
            setupDir = fullfile(folder, [obj.fileName '.set']);
            if exist(setupDir, 'file') == 2
                setupText = fileread(setupDir);
            else
                warning('no setup file found for this measurement');
            end
            
            obj.channels = channel.empty(numel(chNames),0);
            % channel read
            for s=1:numel(chNames)
                fileName = chNames{s};
                fullFileName = fullfile(folder, fileName);
                
                % Check the file exists.
                if( ~exist(fullFileName, 'file') )
                    error('The file does not exist: %s', fullFileName);
                end
                
                % Open the file.
                fileID = fopen( fullfile(folder, fileName), 'r');
                if (fileID == -1)
                    error('The file exists, but could not be opened: %s', fullFileName);
                end
                
                % Read the text header into a variable called 'h'. The loop reads the
                % file character-by-character into h until h finishes with the
                % characters ":CURVE #" or ":CURV #".
                h = '';
                while( isempty( regexp(h, ':CURVE? #', 'once') ) )
                    % If the end of the file has been reached something is wrong.
                    if( feof(fileID) )
                        error('The end of the file %s was reached whilst still reading the header. This suggests that it is not a Tektronix ISF file.', fileName);
                    end
                    c = char(fread(fileID, 1) );
                    h = [h, c];
                end
                
                if s==1
                    obj.waveform_type       = char(regexp(h, 'WFMTYP?E?\s+(.*?)\s*[;:]', 'once', 'tokens'));
                    obj.point_format        = char(regexp(h, 'PT_FM?T?\s+(.*?)\s*[;:]', 'once', 'tokens'));
                    obj.horizontal_units    = char(regexp(h, 'XUNI?T?\s+"*(.*?)"*[;:]', 'once', 'tokens'));
                    obj.horizontal_scale    = str2double(regexp(h, 'HSCAL?E?\s+([-\+\d\.eE]+)', 'once', 'tokens'));
                    obj.horizontal_delay    = str2double(regexp(h, 'HDELA?Y?\s+([-\+\d\.eE]+)', 'once', 'tokens'));
                    obj.sample_interval     = str2double(regexp(h, 'XINC?R?\s+([-\+\d\.eE]+)', 'once', 'tokens'));
                    obj.record_length       = 'Unknown (ISF)';
                    obj.gating              = 'Unknown (ISF)';
                    obj.gating_min          = 'Unknown (ISF)';
                    obj.gating_max          = 'Unknown (ISF)';
                    obj.sample_length       = str2double(regexp(h, 'NR_PT?\s+(\d+)', 'once', 'tokens'));
                else
                    if obj.waveform_type    ~= char(regexp(h, 'WFMTYP?E?\s+(.*?)\s*[;:]', 'once', 'tokens'))
                        error('Waveform Type does not match');
                    end
                    if obj.point_format     ~= char(regexp(h, 'PT_FM?T?\s+(.*?)\s*[;:]', 'once', 'tokens'));
                        error('Point Format does not match');
                    end
                    if obj.horizontal_units ~= char(regexp(h, 'XUNI?T?\s+"*(.*?)"*[;:]', 'once', 'tokens'));
                        error('Horizontal Units do not match');
                    end
                    if obj.horizontal_scale ~= str2double(regexp(h, 'HSCAL?E?\s+([-\+\d\.eE]+)', 'once', 'tokens'));
                        error('Horizontal Scale does not match');
                    end
                    if obj.horizontal_delay ~= str2double(regexp(h, 'HDELA?Y?\s+([-\+\d\.eE]+)', 'once', 'tokens'));
                        error('Horizontal Delay does not match');
                    end
                    if obj.sample_length    ~= str2double(regexp(h, 'NR_PT?\s+(\d+)', 'once', 'tokens'));
                        error('Sample Length does not match');
                    end
                    if obj.sample_interval   ~= str2double(regexp(h, 'XINC?R?\s+([-\+\d\.eE]+)', 'once', 'tokens'));
                        error('Sample Interval does not match');
                    end
                end
                
                % In addition, some header fields are described in the Programmer
                % Manual, but do not seem to appear in any of my files: XMULT, XOFF,
                % XZERO, ZMULT, ZOFF, ZUNIT and ZZERO.
                
                % Check that at least some part of the header was parsed.
                if isempty(str2double(regexp(h, 'BYT_NR?\s+(\d+)', 'once', 'tokens')))
                    warn('Failed to read some part of, or possibly all of, the header in the file %s.', fileName);
                end
                % Calculate the horizontal (x) and vertical (y) values. These equations
                % are given on page 2-171 of the Programmer Manual.
                n = (1:obj.sample_length)';
                if s==1
                    obj.time = (obj.sample_interval * (n - str2double(regexp(h, 'PT_OF?F?\s+([-\+\d\.eE]+)', 'once', 'tokens'))))' - obj.horizontal_delay;
                end
                % read channel signals
                obj.channels(s) = channel.isfreadSignal(obj,fileName,fileID,h);
                
                % read setup file
                if exist('setupText','var')
                    obj.channels(s) = obj.channels(s).isfreadSetup(setupText);
                end
                
                % Close the file
                fclose(fileID);
            end
            
            % -------- math read -----------
            for s=1:numel(mthNames)
                fileName = mthNames{s};
                fullFileName = fullfile(folder, fileName);
                
                % Check the file exists.
                if( ~exist(fullFileName, 'file') )
                    error('The file does not exist: %s', fullFileName);
                end
                
                % Open the file.
                fileID = fopen( fullfile(folder, fileName), 'r');
                if (fileID == -1)
                    error('The file exists, but could not be opened: %s', fullFileName);
                end
                
                % Read the text header into a variable called 'h'. The loop reads the
                % file character-by-character into h until h finishes with the
                % characters ":CURVE #" or ":CURV #".
                h = '';
                while( isempty( regexp(h, ':CURVE? #', 'once') ) )
                    % If the end of the file has been reached something is wrong.
                    if( feof(fileID) )
                        error('The end of the file %s was reached whilst still reading the header. This suggests that it is not a Tektronix ISF file.', fileName);
                    end
                    c = char(fread(fileID, 1) );
                    h = [h, c];
                end
                % create dynamic math property
                if ~isprop(obj,'math')
                    obj.addprop('math');
                end
                
                
                obj.math(s).waveform_type       = char(regexp(h, 'WFMTYP?E?\s+(.*?)\s*[;:]', 'once', 'tokens'));
                obj.math(s).point_format        = char(regexp(h, 'PT_FM?T?\s+(.*?)\s*[;:]', 'once', 'tokens'));
                obj.math(s).horizontal_units    = char(regexp(h, 'XUNI?T?\s+"*(.*?)"*[;:]', 'once', 'tokens'));
                obj.math(s).horizontal_scale    = str2double(regexp(h, 'HSCAL?E?\s+([-\+\d\.eE]+)', 'once', 'tokens'));
                obj.math(s).horizontal_delay    = str2double(regexp(h, 'HDELA?Y?\s+([-\+\d\.eE]+)', 'once', 'tokens'));
                obj.math(s).sample_interval     = str2double(regexp(h, 'XINC?R?\s+([-\+\d\.eE]+)', 'once', 'tokens'));
                obj.math(s).record_length       = 'Unknown (ISF)';
                obj.math(s).gating              = 'Unknown (ISF)';
                obj.math(s).gating_min          = 'Unknown (ISF)';
                obj.math(s).gating_max          = 'Unknown (ISF)';
                obj.math(s).sample_length       = str2double(regexp(h, 'NR_PT?\s+(\d+)', 'once', 'tokens'));
                
                
                % In addition, some header fields are described in the Programmer
                % Manual, but do not seem to appear in any of my files: XMULT, XOFF,
                % XZERO, ZMULT, ZOFF, ZUNIT and ZZERO.
                
                % Check that at least some part of the header was parsed.
                if isempty(str2double(regexp(h, 'BYT_NR?\s+(\d+)', 'once', 'tokens')))
                    warn('Failed to read some part of, or possibly all of, the header in the file %s.', fileName);
                end
                
                if s==1
                    obj.math(s).time = (obj.math(s).sample_interval * (n - str2double(regexp(h, 'PT_OF?F?\s+([-\+\d\.eE]+)', 'once', 'tokens'))))';%  - obj.math(s).horizontal_delay;
                end
                % read channel signals
                obj.math(s).channel = channel.isfreadSignal(obj,fileName,fileID,h);
                % Close the file
                fclose(fileID);
                %
            end
        end
        
        function obj = wfmread(file, verbose,varargin)
            if(~exist('verbose','var'));verbose=-1;warn('All underlying functions are executed in verbose mode');end;
            if ~exist('file', 'var')
                error('No file name, directory or pattern was specified.');
            end
            sizeRowVarargin = cellfun(@(c) size( c,1), varargin, 'UniformOutput', false);
            appIndex= find(cellfun(@(x)~isempty(strfind("app",x)), varargin( [sizeRowVarargin{:}]==1)));
            % Check whether file is a folder.
            if(isa(file, 'char' ) && exist(file, 'dir') )
                folder = file;
                % Get a list of all files that have the extension '.wfm' or '.WFM'.
                files = [ dir(fullfile(folder, '*.wfm')) ];
                
            elseif appIndex==0
                % The pattern is not a folder, so must be a file name or a pattern with
                % wildcards (such as 'Traces/TEK0*.ISF').
                [folder, ~, ~] = fileparts(file);
                % Get a list of all files and folders which match the pattern...
                filesAndFolders = dir(file);
                % ...then exclude the folders, to get just a list of files.
                files = filesAndFolders(~[filesAndFolders.isdir]);
            end
            if isempty(appIndex) || appIndex==0
                fileNames = {files.name};
                datetimes = datestr([files.datenum]);
            end
            obj = scope('Unknown (WFM)');
            obj.firmware_version = 'Unknown (WFM)';
            
            if appIndex>0
                folder = varargin{appIndex +1};
                
                
                chNames = file(~cellfun('isempty',(regexpi(file,'CH'))));
                mthNames = file(~cellfun('isempty',(regexpi(file,'Math'))));
                
                obj.fileName = extractBefore(file{1},'_Ch');
            elseif numel(fileNames)==0
                error('The pattern did not match any file or files: %s', file);
            elseif numel(fileNames) >1 &&  contains(input('Extract one measurement: ','s'),["y","Y","yes","j","ja"])
                fileInDir(:,2) = fileNames';
                fileInDir(:,1)= num2cell(1:numel(fileNames))';
                table(fileInDir)
                fileNr =input('Measurement ID: ');
                fName = extractBefore(fileInDir(fileNr,2),'_Ch');
                fileIDs=~cellfun('isempty',regexp(fileNames,fName));
                fileNames= fileNames(fileIDs);
                
                chNames = fileNames(~cellfun('isempty',(regexpi(fileNames,'CH'))));
                mthNames = fileNames(~cellfun('isempty',(regexpi(fileNames,'MaTH'))));
                
                obj.fileName = fName{1};
            end
            
            obj.channels = channel.empty(numel(chNames),0);
            for s=1:numel(chNames)
                fileName = chNames{s};
                fullFileName = fullfile(folder, fileName);
                
                % Check the file exists.
                if( ~exist(fullFileName, 'file') )
                    error('The file does not exist: %s', fullFileName);
                end
                
                % Read the file.
                [y, t, info] = wfm2read(fullFileName);
                if s==1
                    obj.waveform_type       = info.versioning_number;
                    obj.sample_interval     = 1/info.samplingrate;
                    obj.horizontal_units =  deblank(extractBefore(info.tunit,'!'));
                    obj.record_length       = 'Unknown (WFM)';
                    obj.gating              = 'Unknown (WFM)';
                    obj.gating_min          = 'Unknown (WFM)';
                    obj.gating_max          = 'Unknown (WFM)';
                    obj.sample_length       = info.nop;
                    obj.time = t';
                end
                
                % create signal
                obj.channels(s) = channel.wfmreadSignal(fileName,y,info);
            end
            
            % -------- Math signals ---------
            for s=1:numel(mthNames)
                fileName = mthNames{s};
                fullFileName = fullfile(folder, fileName);
                
                % Check the file exists.
                if( ~exist(fullFileName, 'file') )
                    error('The file does not exist: %s', fullFileName);
                end
                
                % Read the file.
                [y, t, info] = wfm2read(fullFileName);
                
                if ~isprop(obj,'math')
                    obj.addprop('math');
                end
                
                obj.math(s).waveform_type       = info.versioning_number;
                obj.math(s).sample_interval     = 1/info.samplingrate;
                obj.math(s).horizontal_units =  deblank(extractBefore(info.tunit,'!'));
                obj.math(s).record_length       = 'Unknown (WFM)';
                obj.math(s).gating              = 'Unknown (WFM)';
                obj.math(s).gating_min          = 'Unknown (WFM)';
                obj.math(s).gating_max          = 'Unknown (WFM)';
                obj.math(s).sample_length       = info.nop;
                obj.math(s).time = t';
                
                
                % create signal
                obj.math(s).channel = channel.wfmreadSignal(fileName,y,info);
                if strcmp(obj.math(s).channel.name,'math')
                    obj.math(s).channel.name =[obj.math(s).channel.name + char(s)];
                end
                
            end
        end
        
        function obj = csvread(file,channels,verbose,retime)
            if(~exist('verbose','var'));verbose=-1;warn('All underlying functions are executed in verbose mode');end;
            if(~exist('retime','var'));retime=1;end;
            if(~exist('channels','var'));channels={'CH1','CH2','CH3','CH4'};end;
            if(isempty(channels));channels={'CH1','CH2','CH3','CH4'};end;
            if(~exist('file','var'))
                % Ask user for file name
                % check if file exists
            end
            fid = fopen(file);
            cell = textscan(fid, '%*s %s %*s %f', 1, 'delimiter', {',','\n'}, 'headerlines', 0);
            obj = scope(cell{1}{1});
            if strcmp(obj.model,'DPO4054B')
                %% DPO4054B
                obj.firmware_version = cell{2};
                if(verbose);disp('Reading waveform of Tektronix DPO4054B');end;
                DPO4054B_supp_fw = [3.16 3.18 3.20];
                if sum(obj.firmware_version==DPO4054B_supp_fw) == 0
                    if(verbose);warn(['Unsupported version (' num2str(obj.firmware_version) ')']);end;
                    fclose(fid);
                    return;
                end
                cell = textscan(fid, '%s', 8, 'delimiter', {'\n'}, 'headerlines', 1);cell=cell{1};
                [ans cell] = strtok(cell,',');
                cell = strtok(cell,',');
                obj.waveform_type       = cell{1};
                obj.point_format        = cell{2};
                obj.horizontal_units    = cell{3};
                obj.horizontal_scale    = str2double(cell{4});
                obj.horizontal_delay    = str2double(cell{5});
                obj.sample_interval     = str2double(cell{6});
                obj.record_length       = str2double(cell{7});
                obj.gating              = cell{8};
                
                cell = textscan(obj.gating , '%f %f', -1, 'delimiter', {' to ','%'});cell=cell{1};
                obj.gating_min          = cell(1);
                obj.gating_max          = cell(2);
                
                
                obj.sample_length = obj.record_length*((obj.gating_max-obj.gating_min)/100);
                
                
                cell = textscan(fid, '%s', 10, 'delimiter', {'\n'}, 'headerlines', 0);cell=cell{1};
                cell_channels = textscan(cell{10}, '%s', -1, 'delimiter', {','});cell_channels=cell_channels{1};
                fclose(fid);
                data = csvread(file,21,0)';
                obj.time = data(1,1:obj.sample_length);
                if(retime)
                    %retime
                    [minV,minID] = min(abs(obj.time));
                    low = minV-(minID-1)*obj.sample_interval;
                    high = low + obj.sample_length*obj.sample_interval;
                    obj.time = linspace(low,high,obj.sample_length);
                end
                
                k=1;
                for j=1:length(channels)
                    for i=1:length(cell_channels)
                        if(strcmp(channels{j},cell_channels{i}))
                            obj.channels{k}.name = cell_channels{i};
                            obj.channels{k}.value = data(i,1:obj.sample_length);
                            k=k+1;
                        end
                    end
                end
            elseif strcmp(obj.model,'DPO2024')
                %% DPO2024
                obj.firmware_version = cell{2};
                if(verbose);disp('Reading waveform of Tektronix DPO2024');end;
                if obj.firmware_version ~= 1.52
                    if(verbose);warn(['Unsupported version (' num2str(obj.firmware_version) ')']);end;
                    fclose(fid);
                    return;
                end
                cell = textscan(fid, '%*s %s %*s', 7, 'delimiter', {',','\n'}, 'headerlines', 1);cell=cell{1};
                obj.point_format        = cell{1};
                obj.horizontal_units    = str2double(cell{2});
                obj.horizontal_scale    = str2double(cell{3});
                obj.sample_interval     = str2double(cell{4});
                obj.filter_frequency    = str2double(cell{5});
                obj.record_length       = str2double(cell{6});
                obj.gating              = cell{7};
                
                cell = textscan(obj.gating , '%f %f', -1, 'delimiter', {' to ','%'});cell=cell{1};
                obj.gating_min          = cell(1);
                obj.gating_max          = cell(2);
                
                
                obj.sample_length = obj.record_length*((obj.gating_max-obj.gating_min)/100);
                
                
                cell = textscan(fid, '%s', 6, 'delimiter', {'\n'}, 'headerlines', 0);cell=cell{1};
                cell_channels = textscan(cell{6}, '%s', -1, 'delimiter', {','});cell_channels=cell_channels{1};
                fclose(fid);
                data = csvread(file,16,0)';
                obj.time = data(1,1:obj.sample_length);
                if(retime)
                    %retime
                    [minV,minID] = min(abs(obj.time));
                    low = minV-(minID-1)*obj.sample_interval;
                    high = low + obj.sample_length*obj.sample_interval;
                    obj.time = linspace(low,high,obj.sample_length);
                end
                
                k=1;
                for j=1:length(channels)
                    for i=1:length(cell_channels)
                        if(strcmp(channels{j},cell_channels{i}))
                            obj.channels{k}.name = cell_channels{i};
                            obj.channels{k}.value = data(i,1:obj.sample_length);
                            k=k+1;
                        end
                    end
                end
            else
                if(verbose);warn(['Unsupported scope (' num2str(obj.model) ')']);end;
                fclose(fid);
            end
        end
    end
end

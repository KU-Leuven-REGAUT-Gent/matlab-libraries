classdef channel < dynamicprops & matlab.mixin.Copyable
    properties
        nr,
        name, ...
%             filter_frequency, ...
%             record_length, ...
%             sample_length, ...
%             coupling, ...
%             probe_attenuation, ...
            vertical_unit, ...
%             vertical_offset, ...
%             vertical_scale, ...
%             vertical_position, ...
            value, ...
            gain, ...
            forcedRange, ...
            propagationDelay, ...
            invert, ...
            deskew, ...
            bandwidth, ...
            coupling, ...
            offset, ...
            position, ...
            scale, ...
            termination, ...
            label      
    end
    
    methods
        function obj = channel(nr,name)
            obj.nr = str2double(nr);
            obj.name = name;
        end
        
        function obj = decodeChannelPN(obj,s,i,verbose)
            if   ~isprop(obj,'pn')
                obj.addprop('pn');
            end
            obj.pn = eth.empty(1,0);
           scopeTemp =  copy(s);
%            scopeTemp = scopeTemp.copy(s);
            obj.pn = eth.scoperead(scopeTemp,i,verbose);
        end
        
        function [freq_axis, fft_result,N] = advancedFFT(obj,scopeObj,scale,window,gatePosition,gateDuration,verbose)
            
            recordStart = gatePosition-gateDuration/2;
            recordEnd = gatePosition+gateDuration/2;
            periodExtracted = zeros(1,length(scopeObj.time));
            t1 = find(scopeObj.time>=recordStart,1,'first');
            t2= find(scopeObj.time>=(recordEnd-1e-16),1,'first');
            periodExtracted(t1:t2) = 1;
            if exist('verbose','var') && verbose
                figure
                hold on
                plot(scopeObj.time,obj.value)
                plot(scopeObj.time,periodExtracted);
                hold off
            end
            
            sizeData = length(obj.value(t1:t2-1));
            N=sizeData;%2^(nextpow2( sizeData));
            Fs=1/scopeObj.sample_interval;
            
            if sizeData< N
                data =zeros(1,N);
                data(1:sizeData)= obj.value(t1:t2-1);
            else
                data= obj.value(t1:t2-1);
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
            Fs=1/scopeObj(1).sample_interval;
            spectr_res= Fs/N;
            max_freq_possible = Fs/2;
            freq_axis = (Fs*(0:(N/2))/N)';
            
            fft_res = fft(data,N);
            fft_abs= (abs(fft_res/N));
            fft_result=fft_abs(1:N/2+1);
            fft_result(2:end-1) = 2*fft_result(2:end-1);
            if strcmp(scale, 'db')
                fft_result = mag2db(fft_result);
            end
            
        end
        
        function [freq_axis, fft_result,N] = fft(obj,sample_interval,scale,window)
            N= 2^nextpow2( length(obj.value));
            data = zeros(1,N);
            data(1:length(obj.value)) = obj.value;
            switch window
                case 'hanning'
                    %Characteristics:
                    % -     Better frequency,
                    % -     Poorer magnitude accuracy than Rectangular.
                    % -     Hanning has slightly poorer frequency resolution than Hamming.
                    % Best for:
                    % -     Sine, periodic, and narrow-band random noise.
                    % -     Transients or bursts where the signal levels before and after the event are significantly different
                    fftWindow = hanning(N);
                case 'hamming'
                    %Characteristics:
                    % -     Better frequency, ,
                    % -     poorer magnitude accuracy than Rectangular
                    % -     THamming has slightly better frequency resolution than Hanning
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
            Fs=1/sample_interval;
            spectr_res= Fs/N;
            max_freq_possible = Fs/2;
            freq_axis = (Fs*(0:(N/2))/N)';
            
            fft_res = fft(data,N);
            fft_abs= (abs(fft_res/N));
            fft_result=fft_abs(1:N/2+1);
            fft_result(2:end-1) = 2*fft_result(2:end-1);
            if strcmp(scale, 'db')
                fft_result = mag2db(fft_result);
            end
            
        end
        
        function obj = isfreadSetup(obj,setupText)
            
           
            % Get parameter values
            obj.gain = obj.getSetupValues(setupText,'PROBE:GAIN');
            obj.forcedRange = obj.getSetupValues(setupText,'PROBE:FORCEDRANGE');
            obj.propagationDelay = obj.getSetupValues(setupText,'PROBE:PROPDELAY'); % propagation delay of probe (doesn't adjust the data). 
            obj.invert = obj.getSetupValues(setupText,'INVERT');
            obj.deskew = obj.getSetupValues(setupText,'DESKEW'); % Time that signal is shifted in time.
            obj.bandwidth = obj.getSetupValues(setupText,'BANDWIDTH');
            obj.coupling = obj.getSetupParameterString(setupText,'COUPLING');
            obj.offset = obj.getSetupValues(setupText,'OFFSET');
            obj.position = obj.getSetupValues(setupText,'POSITION');
            obj.scale = obj.getSetupValues(setupText,'SCALE');
            obj.termination = obj.getSetupValues(setupText,'TERMINATION');
            obj.label = obj.getSetupParameterString(setupText,'LABEL');     
            
            % Processing
            if obj.deskew ~= 0
                warn( [obj.name ': A deskew of ' num2str(obj.deskew) ' s is configured inside the scope. The signal in MATLAB is already shifted with this time'])
            end          
        end
        function plot(obj,varargin)
            [time,xLimits, title] = obj.splitVarargin(varargin);
              % declaration ylabel
                        switch obj.vertical_unit
                            case 'V'
                                yText = ["Voltage [" + obj.vertical_unit  + "]"] ;
                            case 'A'
                                yText = ["Current [" + obj.vertical_unit  + "]"] ;
                        end
                        if ~isempty(xLimits)

                            startID = find(time >= xLimits(1),1,'first');
                            xLimits(1) = time(startID);
                            endID = find(time <= xLimits(2),1,'last');
                            xLimits(2) = time(endID);
                        else
                            
                        startID = 1;
                        endID = numel(time);
                        end
                        
                        
                        plot(time(startID:endID),obj.value(startID:endID),'LineWidth',2,'Color',pltColor(obj.nr))
                        
                        ylabel(yText)
                        
                        if title ~= ""
                            
                            titlePlot  = obj.name + " - " +  title;
                        else
                            titlePlot = obj.name;
                        end
                        try
                            subtitle(titlePlot);
                        catch
                            suptitle(tiltePlot);
                        end
                        
                        if ~isempty(xLimits)
                            xlim(xLimits);
                        else
                            xlim([time(1), time(end)]);
                        end
        end
    end
    
    methods (Access = private)
        function value = getSetupValues(obj,setupText,parameter)
            regChannel = ['\w*:' upper(obj.name) ':\w*'];
            parameterString = regexp(setupText,[regChannel parameter '[^\n\r]+'],'match');
            value = str2double(extractAfter(parameterString,' '));
        end
        function value = getSetupParameterString(obj,setupText,parameter)
            regChannel = ['\w*:' upper(obj.name) ':\w*'];
            parameterString = regexp(setupText,[regChannel parameter '[^\n\r]+'],'match');
            value = string(extractAfter(parameterString,' '));
        end
        
    end
  
    methods (Static)
      
         function [time, xLimits, titles, savePlot, saveName] = splitVarargin(varargin)          
            varargin = varargin{1};
            if(numel(varargin) >= 1)
                while ~isempty(varargin)
                    if isempty(varargin{1})
                        varargin(1) = [];
                    elseif(ischar(varargin{1}))
                        switch lower(varargin{1})
                            case 'limit'
                                xLimits = varargin{2};
                                varargin(1:2) = [];
                            case 'time'
                                time = varargin{2};
                                varargin(1:2) = [];
                            case 'title'
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
        
            
        function obj = isfreadSignal(scopeObj,fileName,fileID,h)
                BYT_N = str2double(regexp(h, 'BYT_NR?\s+"*(.*?)"*[;:]', 'once', 'tokens'));
                BIT_N = str2double(regexp(h, 'BIT_NR?\s+"*(.*?)"*[;:]', 'once', 'tokens'));
                
                % The next few characters in the file give the number of bytes in the
                % waveform data. The first digit, referred to as 'x' on page 2-60 of
                % the Programmer Manual, gives the number of bytes that immediately
                % follow giving the value 'y', where 'y' is the number of bytes in the
                % waveform. The manual explains it better than I can.
                xBytes = str2double(char(fread(fileID, 1)));
                yBytes = str2double(char(fread(fileID, xBytes)));
                
                % For some reason there is an offset of 1 byte in reading the data
                % files. I don't know why, but I found I could fix it by moving the
                % file position back by one byte.
                fseek(fileID, -1, 'cof');
                
                % Read the waveform.
                % For some oscilloscopes it may be necessary to add 'ieee-be' to the
                % fread statement below. See the comments here:
                % http://www.mathworks.co.uk/matlabcentral/fileexchange/6247-isfread
                if(BYT_N == 2 && BIT_N == 16)
                    [binaryData, count] = fread(fileID, yBytes/2, 'int16');
                elseif (BYT_N == 1 && BIT_N == 8)
                    [binaryData, count] = fread(fileID, yBytes, 'int8');
                else
                    error(['BYT_N ' num2str(BYT_N) ' BIT_N ' num2str(BIT_N) ' - Unknown ISF structure']);
                end
                % Check that the expected number of points have been read.
                if(count ~= scopeObj.sample_length)
                    error('According to the header, the file %s contains %d points, but only %d were read.', fileName, scopeObj.sample_length, count);
                end
                
                % Check that there is no leftover data. I found that there generally
                % is.
                if( ~feof(fileID) )
                    warn('All expected data was read from %s, but there still appears to be data remaining.', fileName);
                end
                channelName = char(extractBetween(h,'WFID "',','));
                nr =extractAfter(lower(channelName),'ch');
                obj = channel(nr,channelName);
                % Calculate the vertical (y) values. These equations
                % are given on page 2-171 of the Programmer Manual.
                obj.value = (str2double( regexp(h, 'YZER?O?\s+([-\+\d\.eE]+)', 'once', 'tokens')) + str2double(regexp(h, 'YMUL?T?\s+([-\+\d\.eE]+)', 'once', 'tokens')) * (binaryData - str2double(regexp(h, 'YOFF?\s+([-\+\d\.eE]+)', 'once', 'tokens'))))';
                obj.vertical_unit = char(regexp(h, 'YUNI?T?\s+"*(.*?)"*[;:]', 'once', 'tokens'));
        end
        
        function obj = wfmreadSignal(fileName,y,info)
            if contains(lower(fileName),'ch')
                nr =extractBetween(lower(fileName),'ch','.wfm');
                obj = channel(nr{1}, char(strcat('CH' , nr{1})));    
            elseif contains(lower(fileName),'math')
                nr = extractBetween(lower(fileName),'math','.wfm');
                obj = channel( nr{1},char(strcat('Math' ,nr{1}))); 
            end
            
                obj.value = y';
                obj.vertical_unit =strcat(info.yunit');
        end
        
    
    end
end
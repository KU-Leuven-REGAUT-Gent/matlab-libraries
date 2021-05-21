classdef channel < dynamicprops
    properties
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
            value
        
    end
    methods
        function obj = channel(name)
            obj.name = name;
        end
        
        function obj = decodeChannelPN(obj,scope,i)
            if   ~isprop(obj,'pn')
                obj.addprop('pn');
            end
            obj.pn = eth.empty(1,0);
            obj.pn = eth.scoperead(scope,i,1);
        end
end
    methods (Static)
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
                
                
                obj = channel(char(extractBetween(h,'WFID "',',')));
                % Calculate the vertical (y) values. These equations
                % are given on page 2-171 of the Programmer Manual.
                obj.value = (str2double( regexp(h, 'YZER?O?\s+([-\+\d\.eE]+)', 'once', 'tokens')) + str2double(regexp(h, 'YMUL?T?\s+([-\+\d\.eE]+)', 'once', 'tokens')) * (binaryData - str2double(regexp(h, 'YOFF?\s+([-\+\d\.eE]+)', 'once', 'tokens'))))';
                obj.vertical_unit = char(regexp(h, 'YUNI?T?\s+"*(.*?)"*[;:]', 'once', 'tokens'));
        end
        
        function obj = wfmreadSignal(fileName,y,info)
            if contains(fileName,'Ch')
                obj = channel( char(strcat('CH' , extractBetween(fileName,'Ch','.wfm'))));    
            elseif contains(fileName,'Math')
                obj = channel( char(strcat('Math' , extractBetween(fileName,'Math','.wfm')))); 
            end
            
                obj.value = y';
                obj.vertical_unit =strcat(info.yunit');
        end
        
    
    end
end
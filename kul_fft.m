%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%                           SCOPE CLASS                              %%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %                                                                      %
%  %    Author: Dimitri De Schuyter                                         %
%  %    Company: KU Leuven                                                %
%  %    Contact: dimitri.deschuyter@kuleuven.be      %
%  %    Version: 1.0                                                      %
%  %                                                                      %
%  %    FFT that results a frequency axis, spectrum and windows size                        %
%  %    Readable files: isf, csv                                          %
%  %                                                                      %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %        FUNCTIONS (static)                 *Object creation*          %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %                                                                      %
%  %    usage: [freq_axis, fft_result,N] =kul_fft(obj, scale,window.)                 %
%  %                                                                      %
%  %    obj needs obj.value, obj.sample_interval,   %
%  %      scale                    'db' or empty        %
%  %        window                hanning, hamming, blackmanharris, rectangle, gaussian       %

function [freq_axis, fft_result,N] = kul_fft(obj,scale,window)
            N= 2^nextpow2( length(obj.value));
            data = zeros(1,N);
            data(1:length(obj.value)) = obj.value;
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
            freq_axis(1,:) = (Fs*(0:(N/2))/N)';
            
            fft_res = fft(data,N);
            fft_abs= (abs(fft_res/N));
            fft_result(1,:)=fft_abs(1:N/2+1);
            fft_result(1,2:end-1) = 2*fft_result(2:end-1);
             if strcmp(scale, 'db')
                   fft_result = mag2db(fft_result);
             end
            
        end
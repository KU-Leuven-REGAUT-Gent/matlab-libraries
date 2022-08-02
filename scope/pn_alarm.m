close all;
addpath('Libraries');

objEth = eth.pcapread('pn_alarm.pcapng',0);

j = 0;
for i=1:length(objEth)
    try
        if(objEth(i).EtherTypeSpecificData.PNIO_FrameID == 'FE01')
            j = j+1;
            objAlarm(j)=objEth(i);
        end
    catch
        warning('No PNIO_FrameID');
    end
end

       
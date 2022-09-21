classdef pcapMeasurement
    properties
        fileName
        port
    end
    
    methods
        function obj = pcapMeasurement(file)
            % reads the pcap file and creates a pcapMeasurement object that
            % stores the filename and creates a port property that contains
            % an array that corresponds with the interface ids that are
            % inside the pcap. All packets of the specific interface id
            % will be stored in to the corresponding port.
            LastFolderIndex = find(file=='\',1,'last')+1;
            obj.fileName = cell2mat(extractBetween(file,LastFolderIndex,'.pcap'));
            obj.port = eth.pcapread(file,0);
            
        end
        
        function [fig, table] = table(obj, fig)
            % generates a table of all packets in obj.ports.
            % Implemented protocols:
            % -  RT class 2
            % -  Alarm Low
            % -  Alarm High
            % -  MRP
            % -  preempted

            if(~exist('fig','var'));fig = uifigure('Position',[300 100 1200 800]);end; 
            % create the data
            
            
            colorgen = @(text,color) ['<html><div style="color:rgb(',num2str(color(1)),',',num2str(color(2)),',',num2str(color(3)),');">',text,'</div></html>'];
            bgcolorgen = @(text,color,bgcolor) ['<html><div style="width:200px;height:15px;padding:2px;color:rgb(',num2str(color(1)),',',num2str(color(2)),',',num2str(color(3)),');background-color:rgb(',num2str(bgcolor(1)),',',num2str(bgcolor(2)),',',num2str(bgcolor(3)),');">',text,'</div></html>'];
            
            % get information of all packets in the port variable
            i=1;
            for portID=1:length(obj.port)
                for packetID = 1:length(obj.port(portID).packet)                    
                    d{i,1} = num2str(portID);
                    d(i,2:7) = struct2cell(obj.port(portID).packet(packetID).getInfo');             
                    i=i+1;
                end
            end
            d =sortrows(d,2);
            
            % Create the column and row names in cell arrays
            columnname = {'Interface ID','Time','Destination address','Source address','Packet description','Type','Data'};
            columnformat = {'char','numeric','char','char','char','char','char'};
            columnwidth = {50,100,120,120,300,200,200};
            % Create the uitable
            tble = uitable(fig,'Data', d,...
                'ColumnName', columnname,...
                'ColumnFormat', columnformat,...
                'ColumnWidth', columnwidth,...
                'RowName',[],...
                'Position',[50 50 fig.Position(3)-100 fig.Position(4)-100],...
                'BackgroundColor',[1 1 1]);
            tble.RowName = 'numbered';
            pnStyle = uistyle('BackgroundColor','green');
            alarmLowStyle = uistyle('BackgroundColor',[0.8500 0.3250 0.0980]);
            alarmHighStyle = uistyle('BackgroundColor','red');
            mrpStyle = uistyle('BackgroundColor',[0.4940 0.1840 0.5560]);
            preemptionStyle = uistyle('BackgroundColor',[0.3010 0.7450 0.9330]);
            
            % Find corresponding protocol packets
            pnPacketsID =find(contains(d(:,5),'RT class 2'));
            alarmLowPacketsID =find(contains(d(:,5),'Alarm Low'));
            alarmHighPacketsID =find(contains(d(:,5),'Alarm High'));
            mrpPacketsID = find(contains(d(:,5),'MRP'));
            preemptionID = find(contains(d(:,6),'preempt'));
            
            % Adjust the style of each row to the corresponding protocol
            % style
            addStyle(tble,pnStyle,'row', pnPacketsID);
            addStyle(tble,alarmLowStyle,'row', alarmLowPacketsID);
            addStyle(tble,alarmHighStyle,'row', alarmHighPacketsID);
            addStyle(tble,mrpStyle,'row', mrpPacketsID);
            addStyle(tble,preemptionStyle,'row', preemptionID);
        end
        
    end
end


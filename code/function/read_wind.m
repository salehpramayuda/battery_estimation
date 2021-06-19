function [data, speed] = read_wind(filename)
    %% Read wind data from xls file from rp5.ru
    data = readcell(filename);
    data_ = data(8:end,:);

    % check missing data
    miss = 0;
    where = [];
    for i=1:length(data_)
        temp = data_(i,7); temp=temp{1};
        if(ismissing(temp))
            miss = miss+1;
            where = cat(2,where,i);
        end
    end

    % create struct array
    % struct.Date : Date of data
    % struct.Time : Time of data
    % struct.Vector : speed*[east-west, north-south]
    data = repmat(struct('Date','','Time','','Speed',0,'Direction',''),length(data_)-miss,1);
    speed = zeros(1, length(data_)-miss);
    for i = 1:length(data_)

        % extract one data
        j=length(data_)-i+1;
        if(any(j==where))
            j=j-1;
            i=i+1;
        end
        temp = data_(j,1); temp = temp{1};
        [date,time] = strtok(temp,' ');
        data(i).Date = date;
        data(i).Time = strtok(time,' ');
        
        % check wind speed
        spid = data_(j,8); spid = spid{1};
        data(i).Speed = spid;
        speed(i) = spid;
        if(spid==0)
            data(i).Direction = 'north';
        else
            temp = data_(j,7); temp=temp{1};
            while(~strcmp(temp,''))
                [token,temp]=strtok(temp,' ');
            end
            data(i).Direction = token;
        end
    end
    
    % percentage of speed above 9 and 10 m/s
    per_09 = sum(speed>=9)/length(speed)*100;
    per_10 = sum(speed>=10)/length(speed)*100;
    
    % plot wind speed
    figure();
    plot(speed);
    xlabel(['Days', newline, num2str(per_09),'% above or equal 9m/s',...
        newline, num2str(per_10), '% above or equal 10m/s']);
    ylabel('Speed [m/s]');

%         % check wind speed
%         spid = data_(j,8); spid = spid{1};
%         if(spid~=0)
%             % parse wind direction
%             temp = data_(j,7);temp = temp{1};
%             while(~strcmp(temp,''))
%                 [token,temp]=strtok(temp,' ');
%             end
%             [direct1,direct2]=strtok(token,'-');
% 
%             % check direction
%             if(strcmp(direct1,'north')) 
%                 phase = 90/180*pi;
%             elseif(strcmp(direct1,'east')) 
%                 phase = 0;
%             elseif(strcmp(direct1,'south'))
%                 phase = -90/180*pi;
%             else
%                 phase = pi;
%             end
% 
%             % check 2nd direction
% 
%             if(~isempty(direct2))
%                 direct2=strtok(direct2,'-');
%                 if strcmp(direct2,'north')
%                     phase = 1/2*(phase+90/180*pi);
%                 elseif strcmp(direct2,'northeast')
%                     phase = 1/2*(phase+45/180*pi);
%                 elseif strcmp(direct2, 'east')
%                     phase = 1/2*(phase+0);
%                 elseif strcmp(direct2,'southeast')
%                     phase = 1/2*(phase-45/180*pi);
%                 elseif strcmp(direct2, 'south')
%                     phase = 1/2*(phase-90/180*pi);
%                 elseif strcmp(direct2,'southwest')
%                     phase = 1/2*(phase-135/180*pi);
%                 elseif strcmp(direct2, 'west')
%                     phase = 1/2*(phase+pi);
%                 else
%                     phase = 1/2*(phase+135/180*pi);
%                 end
%             end
% 
%             data(i).Vector = spid*[round(cos(phase),10),round(sin(phase),10)];
%             data(i).Speed = spid;
%             
%         else
%             data(i).Speed = spid;
%             data(i).Vector = [0, 0];
%         end

end
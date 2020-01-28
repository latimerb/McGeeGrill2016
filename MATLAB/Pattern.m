function [ svec ] = Pattern(PatternSelect, pulses, pulsew, stimdel, t_steps, dt )
stimtime=zeros(1,max(pulses));
stimtime(1)=stimdel;

%% PATTERNS
% Pattern 1: Constant 33 Hz
if PatternSelect == 1
    for s = 2:max(pulses)
        stimtime(s) = stimtime(s-1)+round((1/33)*(1000/dt)); %time when pulses occur in timesteps
    end
end

% Pattern #2 : DECREASING 
if PatternSelect == 2
    for s=2:max(pulses)
        stimtime(s)=stimtime(s-1)+round((1./(50-((25.4/max(pulses))*s)))*(1000/dt));
    end
end

% Pattern #3 : INCREASING
if PatternSelect ==3 
    for s=2:max(pulses)
        stimtime(s)=stimtime(s-1)+round((1./(24.6+((25.4/max(pulses))*s)))*(1000/dt));
    end
end

% Pattern #4 : RANDOM
if PatternSelect ==4
    for s=2:max(pulses)
        stimtime(s)=stimtime(s-1)+round((1/(33+(randn*7)))*(1000/dt));
    end
end

% Pattern #5 : GAPS 50 %100 ms on 100 ms off
if PatternSelect == 5
    ss=[];
    counter=round(.1/(1/66)); %ms/freq ipi -> 7 on 7 off
    for p=1:2:floor((max(pulses)*2)/counter)
        ss=[ss,counter*(p-1)+1:counter*(p-1)+7];
    end
    
    for s = 2:max(pulses)*2 %66 Hz
        stimtime(s) = stimtime(s-1)+round((1/66)*(1000/dt)); %time when pulses occur in timesteps
    end 
    stimtime=stimtime(ss); 
end  

% Pattern #6 : Doublets 10 (10/50ms)
if PatternSelect ==6
    for s=2:max(pulses)
        if mod(s,2)==0
            %number is even
        stimtime(s) = stimtime(s-1)+(.0102)*(1000/dt); %10 ms IPI
        else
            %pulse number is odd
        stimtime(s) = stimtime(s-1)+(.0501)*(1000/dt); %50 ms IPI
        end
    end
             
end

% Pattern #7: Alternating (20/40ms)
if PatternSelect == 7
    for s=2:max(pulses)
        if mod(s,2)==0
            %number is even
        stimtime(s) = stimtime(s-1)+(.0202)*(1000/dt); %20 ms IPI
        else
            %pulse number is odd
        stimtime(s) = stimtime(s-1)+(.0401)*(1000/dt); %40 ms IPI
        end
    end
end

% Pattern #8: 33 Pauses
if PatternSelect == 8 
        
    del=floor(.2*33);% #pulses to remove
    ss=1:max(pulses);
    for sc=1:(max(pulses)/33)
        ss((sc*33)-del:(sc*33))=0; %delete last 0.2 of every second of pulses
    end
    for s = 2:max(pulses) %66 Hz
        stimtime(s) = stimtime(s-1)+round((1/33)*(1000/dt)); %time when pulses occur in timesteps
    end 
    stimtime=stimtime(ss~=0); %get rid of pulses in pauses
end

% Pattern #9 : 33 Hz Doublets
if PatternSelect == 9
    for s=2:max(pulses)*2
        if mod(s,2)==0
            %number is even
        stimtime(s) = stimtime(s-1)+(.01)*(1000/dt); %10 ms IPI
        else
            %pulse number is odd
        stimtime(s) = stimtime(s-1)+(.02)*(1000/dt); %20 ms IPI
        end
    end
end


%% AFTER Pattern and stimtimes are decided, create stim vector with 1s in
% timesteps
svec=zeros(1,t_steps);
for a = 1:length(stimtime)
svec(stimtime(a)-pulsew:stimtime(a)+pulsew) = 1;
end

end


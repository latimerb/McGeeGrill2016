close all; clear all; clc;

%% ---GENERAL LIAF PROPERTIES------------------------------------%
v_rest = -70;      % mV Resting membrane voltage
v_thresh = -50;    % mV Threshold voltage -> spike
v_spike = 20;      % mV Voltage Vm goes to after crossing Vthresh
v_th = v_thresh;
v = zeros(nn,t_steps);
v(:,1) = v_rest;   % initial condition (resting potential)
u = zeros(nn,t_steps);
v_rest = -65*ones(nn,1); % resting potential (mV) 
tau_m = 10*ones(nn,1);  % membrane time constant (ms) 
R_m = 10*ones(nn,1);  % membrane resistance (Mohm)  
ga = 0.1* ones(nn,1);
gad = zeros(nn,1);

%---ACTION POTENTIAL PARAMETERS---------------------------------%
v_ap = 60; % max membrane potential of AP (mV)
t_ap = 0;  % timestamp of AP (ms)
dur_ap = 0;  % duration of AP (ms)
dur_rp = 1;  % duration of refractory period (ms)
no_ap = 0;  % number of AP
fire = 0;  % boolean (0 = not firing AP, 1 = firing AP)

G = zeros(nn,nn);
Z = zeros(nn,nn);
g = zeros(nn^2,t_steps); % recording values in G
z = zeros(nn^2,t_steps); % recording values in Z
fire_ap = (v(:,1) > v_th);
time_ap = zeros(nn,1)-2*dur_rp;
v_int = zeros(nn,t_steps); 
windowlength= 1000/dt; %1s window

%% ---Time Marching---------------------------------------------%
v_th = v_thresh*ones(nn,t_steps);
for k = 2:t_steps
    % AFFERENT - Generate Pelvic Firing Rate from Bladder Pressure at last timestep
        % Pelvic Afferents: Low Threshold type neuron
        x=PBladder(k-1); %cmH2O
        FRlow = -3E-08*x^5 + 1E-05*x^4 - 0.0015*x^3 + 0.0799*x^2 - 0.5924*x;
           FRlow(FRlow<0)=0;
        PelAffIPI=round(1000*(1/FRlow/dt)); % msec/dt= in timesteps
   
    t_from_ap = t(k) - time_ap;
    ck_rr = (t_from_ap >=0) & (t_from_ap < dur_rp);    
    v_int = v(:,k-1);
    v_int(v_int>v_ap)=v_ap; 
    if(any(fire_ap))
        v_int(fire_ap) = v_rest(fire_ap);
    end       
     ga(:,k) = ga(:,k-1)+(-ga(:,k-1)/tauad)*dt;
     gad = ga(:,k);
  
     % Set Pelvec spiketimes based on Bladder Pressure
     Lastpel=find(Pelvec==1, 1, 'last' ); %last place holder of firing time of pelvec %or max(find(Pelvec==1));
     Pelvec(k)=(k-Lastpel>=PelAffIPI); %set=1 if time elapsed >= to Pel Afferent IPI
  
      % Set PMC firing rate based on window of Pel Afferents from last timestep
      % PMCRate %Hz from previous timestep of preceding window
      PMCIPI= round(1000*(1/PMCRate(k-1)/dt)); 
      LastPMC=find(PMCvec==1,1,'last');
      PMCvec(k)=(k-LastPMC>=PMCIPI); %set = 1 if time elapsed is >= IPI for FR of PMC
     
    i_stim = zeros(nn,1); % reset stim, so it doesn't increase size of vec.
    i_stim(1,stimvec(round(t(k)*10))==1)=i_mag*2; % if stimvec=1, then set i_mag
    i_stim(1,Pudvec(round(t(k)*10))==1)=i_mag*5; % Set Pud firing times
    i_stim(2,Pelvec(k)==1)=i_mag*8; %Set Pel firing times 
       % if k<windowlength % ie before closing the loop happens
        i_stim(2,Pelveca(k)==1)=i_mag*10;
       % end
    i_stim(8,PMCvec(k)==1)=i_mag*8;
    
    % Calculate Synatpic Current------------------------------------%
    Vsyn = (v_int*ones(1,nn))';
    Usyn = (u(:,k-1)*ones(1,nn));
    isyn = (Gbar.*G.*(Vsyn-Erev))'*ones(nn,1);

    % CALCULATE neuron membrane potential---------------------------% 
    dvdt = (1./tau_m).*( -(v_int - v_rest).*(1+(R_m.*gad)) +R_m.*(i_stim - isyn)) ;
         
    dGdt = -G./Tau_d + Z;
    dZdt = -Z./Tau_r + gsyn.*Usyn;
     
    v(:,k) = v_int + dvdt*dt;
    v(ck_rr,k) = v_rest(ck_rr);
    G = G + dGdt*dt;
    Z = Z + dZdt*dt;
    g(:,k) = G(:);
    z(:,k) = Z(:);
    % When action potential is fired:
    fire_ap = (v(:,k) > v_th(:,k));
    time_ap(fire_ap) = t(k);
    v(fire_ap,k) = v_ap;
    u(fire_ap,k) = 1/dt;
    ga(fire_ap,k) = ga(fire_ap,k)+gainc;
    
    BladderVol=BladderVols(k);
   % Pelvic Efferent Firing Analysis ----------------------------%
   % Windowing] -> P_Bladder % Average Output firing rate
    if k>windowlength
    % when time > window length, then start updating the frequencies/pressures
    outfire=zeros(1,windowlength+1); % 500ms/dt=5000 time steps = 0.5 sec window
    outfire=outfire(v(7,k-windowlength:k)>=0);  
    OUTFIRE= length(outfire)/(windowlength/(1000/dt)); % Hz
    VolCalc= ((1.5*BladderVol-10)>0)*((1.5*BladderVol-10)); % calculate contribution to Pbl from just volume filling alone
    PBladder(k) = VolCalc +(0.0002*OUTFIRE^3 - 0.033*OUTFIRE^2 + 1.8045*OUTFIRE - 0.503);%equation to convert freq to bladder pressure
    
    % PelAfferents FR to set PMC output rate
    pelfiring=find(v(2,k-windowlength:k)>0); % indices where spike occurs
    Pelipis=pelfiring(2:length(pelfiring))-pelfiring(1:length(pelfiring)-1); % IPI of pel aff in timesteps
    PelvicFR(k)=mean(1./(Pelipis/10000));
    PelPMCsig=(sum(Pelipis>1500)<1); % binary signal to fire PMC
    PMCRate(k)=PelPMCsig*15*(BladderVols(k)>10); % if on, the fire at 20 Hz, if no, at 0 
    
    end
    if mod(k,100)==0
        fprintf('%f percent done\n',100*k/t_steps)
    end
end
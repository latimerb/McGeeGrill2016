close all;clear all;clc;

%% Sacral Spinal Network Model %%
%% Version 11
% McGee and Grill 2016, J Comput Neurosci 

% This code is used to model the simplified neural network model of the
% pudendo-vesical reflex using linear integrate and fire (LIAF) neurons.
% The supraspinal conribution of PAG/PMC are represented by a simple
% transfer function producing on/off response
% Use the following section to specify easy-access inputs.

% Neural Network Structure:
% Pud = Pudendal Afferent Input
% Pel = Pelvic Afferent Input
% IN_d = Dorsal Interneuron
% IN_mn = Medial Interneuron (inhibitory)
% IN_mp = Medial Interneuron (excitatory)
% FB = Feedback Neuron
% Out = Pelvic Output
% PMC = Pontine Micturition Center output (descending on/off signal)

%% Regular and/or Temporal Patterns: loop for both or select 0/1.    
for sss=1
SELECT=1; % 0=Temporal Patterns, 1=Regular Frequency Stimulation
% Turn on output plots of bladder pressure and firing activity
BladderPlot=1; % Plot bladder if =1, suppress plot if =0
% Turn on bladder filling
Filling=1; % Bladder state: 1=ramp/filling, 0=isovolumetric
% Initial bladder volume:
BladderVol= 10; % avg 100%= 15mL, low pressure until 8-10 ml 
%------------------------------------------------------------%
% Time variables:
TMAX=50; % sec
StimLength=5; % sec  
StimDelay=TMAX; % sec do not make less than 2 sec = calc baseline
InitDelay=2; % Initialization delay to account for windowing
%------------------------------------------------------------%
% Multiple Stimulation Testing (i.e., loop through code multiple times to compare
% different frequencies or patterns of stimulation)
if SELECT==1
% WHAT FREQUENCIES DO YOU WANT TO STIMULATE?
FreqCases=[33]; % regular frequency stimulation
FC=length(FreqCases); 
PatternSelect=1;
ModelOut=zeros(FC,14);
Identifier='freqs';
    else % Temporal Patterning
% WHICH PATTERNS DO YOU WANT TO TEST?
PatternSelect=[1];
FC=length(PatternSelect);
ModelOut=zeros(FC,13);
Identifier='patterntest';
end
%------------------------------------------------------------%
%% Network Topology %%
% Run model & analysis for EACH stimulation frequency or pattern (total #loops =FC)
for q=1:FC % q is a counter of which one you are on
% Network Synapse Properties---------------------------------%
Aa=0.45; % a Pel excites IN_d 
Bb=0.6;  % b Pud excites IN_d  
Cc=0.44; % c Pud excites IN_mp
Dd=0.7;  % d Pud excites IN_mn
Ee=.65;  % e IN_mn inhibits Out
Ff=0.6;  % f IN_mp excites Out
Gg=0.8;  % g IN_d excites Out
Hh=1.0;  % h Out excites FB
Ii=0.6;  % i FB inhibits IN_d
Jj=0.33; % j PMC excites Out
%------------------------------------------------------------%
%% Variable Initialization %% 
%---TIME-----------------------------------------------------%
t_max = TMAX*1000; % ms (not steps) 
dt = 0.1;          % ms
t = 0:dt:t_max;    % ms
t_steps = length(t);

%---DEFINE NEURONS--------------------------------------------%
nn_names=['Pud  '; 'Pel  '; 'IN_d '; 'IN_mn'; 'IN_mp'; 'FB   '; 'Out  ';'PMC  '];
% Neurons are arranged in a vector format, each neuron always refers to the
% same index of whatever variable/parameter is needed.
%For example, to reference Neuron #2's name: nn_names(2,:) = Pel
nn = length(nn_names); % total number of neurons = 8 (pudendal has multiple inputs onto interneurons)

%% ---NEURON FIRING ACTIVITY-----------------------------------%
i_dur = t_max;   % duration of stimulus pulse (ms)
i_delay = 5;     % delay before stimulating (ms)
i_mag = 30;      % magnitude of stimulus pusle (nA) %3=40Hz
stimvec = zeros(1,t_steps);
stimdel = StimDelay*10000;   % steps, stim. start delay [1=Pud,2=Pel] 
stimend = stimdel+(StimLength*10000);   %when stim stops %t_steps-10000 if want 1sec before
pulsew = 3;      % timesteps, half pulsewidth 

%---SET PN STIMULATION-----------------------------------------%
% Create vectors of firing times PN Stimulation based on freq/pattern
if SELECT==1 % REGULAR FREQUENCY STIMULATION
    PNSTIM=FreqCases(q);
    stimrate = PNSTIM;   % Hz, Frequency of PN/Pel Stimulation
    stimipi = round((1./stimrate*10000)); %in timesteps
    stimtime=zeros((round((stimend-min(stimdel))/min(stimipi))));
    stimtime(:,1)=stimdel; 
    for s = 2:(round((stimend-min(stimdel))/stimipi))
        stimtime(s) = stimtime(s-1)+stimipi;
        s = s+1;
    end
    for s = 1:length(stimtime)
        stimvec(stimtime(s)-pulsew:stimtime(s)+pulsew) = 1;
    end
else if SELECT==2
        stimrate1=Frequencies(1);
        stimipi1 = round((1./stimrate1*10000)); %in timesteps
        stimtime=zeros((round((stimend1-min(stimdel1))/min(stimipi1))));
        stimtime(:,1)=stimdel; 
        for s = 2:(round((stimend1-min(stimdel1))/stimipi1))
            stimtime(s) = stimtime(s-1)+stimipi1;
            s = s+1;
        end
        for s = 1:length(stimtime)
        stimvec(stimtime(s)-pulsew:stimtime(s)+pulsew) = 1;
        end
        stimrate2=Frequencies(2);
         stimipi2 = round((1./stimrate2*10000)); % in timesteps
        stimtime(:,1)=stimdel2; 
        for s = 2:(round((stimend2-min(stimdel2))/stimipi2))
        stimtime(s) = stimtime(s-1)+stimipi2;
        s = s+1;
        end
        for s = 1:length(stimtime)
        stimvec(stimtime(s)-pulsew:stimtime(s)+pulsew) = 1;
        end
else % TEMPORAL PATTERNS OF STIMULATION *(USES FUNCTION Pattern.m)*
    stimrate = 33;   % Hz, Frequency of PN/Pel 
    PNSTIM= PatternSelect(q);
    stimipi = round((1./stimrate*10000)); % in timesteps
    pulses=1:round((stimend-stimdel)/stimipi);
    stimvec=Pattern(PatternSelect(q), pulses, pulsew, stimdel, t_max, dt);
    stimvec(length(stimvec):t_steps)=0;
    end
end

%---SET PN INPUT FIRING ACTIVITY---------------------------------%
% Define Natural Sensory firing pattern of Pel and Pud
% Pudendal Firing Base Activity
Pudrate=0; % Baseline average of PN activity 
Pudipi= round((1/Pudrate*10000));
Puddel=80;
Pudend=t_steps;
pulsew = 3;      % timesteps, half actual pw
Pudvec = zeros(1,t_steps);
Pudstimtime=zeros((round((Pudend-Puddel)/Pudipi)));
Pudstimtime(1)=Puddel;
for s = 2:(round((Pudend-Puddel)/Pudipi))
    Pudstimtime(s) = Pudstimtime(s-1)+Pudipi+round(randn*500); % stochastic element, 5=st. dev 0.5 ms
    s = s+1;
end
for s = 1:length(Pudstimtime)
Pudvec(Pudstimtime(s)-pulsew:Pudstimtime(s)+pulsew) = 1;
end

%---SET PELVIC INPUT FIRING ACTIVITY---------------------------%
Pelrate=1;  % Initial Pelvic Firing rate before loop takes over 
Pelipi= round((1/Pelrate*10000));
Peldel=25;
Pelend=t_steps;
pulsew = 3;      % timesteps, half actual pw
Pelvec = zeros(1,t_steps);
PMCvec = zeros(1,t_steps);
Pelveca = zeros(1,t_steps);
Pelstimtime=zeros((round((Pelend-Peldel)/Pelipi)));
Pelstimtime(1)=Peldel;
for s = 2:(round((Pelend-Peldel)/Pelipi))
    Pelstimtime(s) = Pelstimtime(s-1)+Pelipi+round(randn*70); % stochastic element, 5=st. dev 0.5 ms
    s = s+1;
end
 for s = 1:length(Pelstimtime)
 Pelveca(Pelstimtime(s)-pulsew:Pelstimtime(s)+pulsew) = 1;
 end
Pelvec(1)=1;
PMCvec(1)=1;
%%---BLADDER FILLING CALCULATION-------------------------------%
PBladder=15*ones(1,t_steps); % cm H2O,needs a few sec to initialize
BladderVols=zeros(1,t_steps);
if Filling==1
% Bladder Volume Infusion --- RAMP ---
BladderVols= (17*t)/t_max; 
else
% Isovolumetric Bladder Volume --- CONSTANT ---
BladderVols= BladderVol*ones(1,t_steps);%     
end

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
%---SET SYNAPTIC PROPERTIES------------------------------------%
Gbar = zeros(nn,nn);  % conductance matrix for neurons
Gbar(2,3)=Aa; % a Pel excites IN_d 
Gbar(1,3)=Bb; % b Pud excites IN_d 
Gbar(1,5)=Cc; % c Pud excites IN_mp 
Gbar(1,4)=Dd; % d Pud excites IN_mn 
Gbar(4,7)=Ee; % e IN_mn inhibits Out
Gbar(5,7)=Ff; % f IN_mp excites Out
Gbar(3,7)=Gg; % g IN_d excites Out
Gbar(7,6)=Hh; % h Out excites FB
Gbar(6,3)=Ii; % i FB inhibits IN_d
Gbar(8,3)=Jj; % j PMC excites Out

Erev = zeros(nn,nn); % default for all excitatory synapses
Erev(4,7) = -80; % inhibitory synapse
Erev(6,3) = -80; % inhibitory synapse

% These modify individual IPSP/EPSP shape
Tau_r = 0.9*ones(nn,nn); % AMPA 
Tau_r(4,7)=1.1; % GABAA
Tau_r(6,3)=1.1; % GABAA
Tau_d = 12.15*ones(nn,nn); % AMPA 
Tau_d(4,7)=10; % GABAA
Tau_d(6,3)=10; % GABAA

gpeak=ones(nn,nn)*0.28; %mS/cm^2 AMPA 
gpeak(4,7)=1.5; % mS/cm^2 GABAA
gpeak(6,3)=1.5; % mS/cm^2 GABAA
 
peaktime = ((Tau_d.*Tau_r)./(Tau_d-Tau_r)).*log(Tau_d./Tau_r);
gsyn = gpeak./((Tau_d.*Tau_r)./(Tau_d-Tau_r)).*(exp(-peaktime./Tau_d)-exp(-peaktime./Tau_r));
tau_th = 40;
gainc = 0.5;
tauad = 35;

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
% set PMC rate to 0 at baseline 
PMCRate=zeros(1,t_steps);

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

%% ---Contraction Analysis ------------------------------------------%
%%---Individual Run Plots -------------------------------------------%
if BladderPlot==1
figure()
hold on;
for ii = 1:nn
    subplot(nn,1,ii),plot(t,v(ii,:),'k');
    ylim([-100,0]);
    xlim([0,10000]);
    if ii==1
    ylabel(char(nn_names(1,:)))
    title([sprintf('%6.2g',Aa,Bb,Cc,Dd,Ee,Ff,Gg,Hh,Ii),'  PN Stim:',sprintf('%6.2g',PNSTIM),' Hz'])
    else if ii==2
    ylabel(char(nn_names(2,:)))
            else if ii==3
            ylabel(char(nn_names(3,:)))
                else if ii==4
                ylabel(char(nn_names(4,:)))
                    else if ii==5
                    ylabel(char(nn_names(5,:)))
                        else if ii==6
                        ylabel(char(nn_names(6,:)))
                            else if ii==7
                            ylabel(char(nn_names(7,:)))
                                else
                                ylabel(char(nn_names(8,:)))
                            end
                        end
                    end
                end
        end
    end
end
end
% PLOT of Afferent inputs, output, and PBladder
figure() 
subplot(6,1,1) % change to 5 to include volume infused
plot(t,v(1,:),'k');
if SELECT==1
    title([sprintf('%6.2g',Aa,Bb,Cc,Dd,Ee,Ff,Gg,Hh,Ii),'  PN Stim:',sprintf('%6.2g',PNSTIM),' Hz'])
else title([sprintf('%6.2g',Aa,Bb,Cc,Dd,Ee,Ff,Gg,Hh,Ii),'      Pattern #',sprintf('%6.3g',PatternSelect(q))])
end
ylabel('Pud Input')
ylim([0 50])
xlim([500 t_max])
subplot(6,1,2)     %
plot(t,v(2,:),'k');
ylabel('Pel Input')
ylim([0 50])
xlim([500 t_max])
subplot(6,1,3)      %
plot(t,v(7,:),'k');
ylim([0 50])
xlim([500 t_max])
ylabel('Pel Output')
subplot(6,1,4)        %
plot(t,smooth(PBladder,1500),'LineWidth',1.5) % Uses Curve Fitting Toolbox
%plot(t,PBladder,'LineWidth',1.5) % Use this instead if toolbox above not available
ylabel('P_B_l_a_d_d_e_r')
ylim([0 50])
xlim([500 t_max])
StimLinex=[stimdel,stimend]; StimLiney=[4,4];
line(StimLinex*dt,StimLiney,'LineWidth', 3,'Color','k')
subplot(6,1,5)
plot(t,BladderVols)
xlim([500 t_max])
ylabel('Volume (mL)')
subplot(6,1,6)
plot(t,v(8,:),'k')
ylabel('PMC')
xlabel('Time (ms)')
end

%% Bladder Pressure and Firing Rate Analysis ---------------------%
%BLADDER PRESSURE
PRE=mean(PBladder(stimdel/2:stimdel));
 STIM=mean(PBladder(stimdel:stimend));
 POST=mean(PBladder(stimend:t_steps));
 ContractionMag=STIM-PRE; 
%FIRING RATE
outfirecont=zeros(1,stimdel+1); %500ms/dt=5000 time steps = 0.5 sec window
    outfirecont=outfirecont(v(7,1:stimdel)>=0);
    OUTFIREconts= length(outfirecont)/(stimdel/(1000/dt)); % Hz
outfirestim=zeros(1,stimend-stimdel+1); %500ms/dt=5000 time steps = 0.5 sec window
    outfirestim=outfirestim(v(7,stimdel:stimend)>=0); %=1?*** CHECK
    OUTFIREstims= length(outfirestim)/((stimend-stimdel)/(1000/dt)); % Hz

q=q+1;
end
end

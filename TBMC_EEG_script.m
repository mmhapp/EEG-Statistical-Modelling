% Use 'Evaluate Selection in Command Window' to run each module of code

% Command to get the EEGLAB window
eeglab redraw

% load the dataset
subjNumber = 4; % Change accordingly, total of 4 subjects
EEG = pop_loadset(['/Users/mikkohapponen/Documents/Tiedostot/Opiskelu/DI/Electroencephalography/EEG Signals/Datasets/Raw/Subject' num2str(subjNumber) 'Raw.set']);


%% Basic EEG handling

% We make this variable for guiding future channel interpolation
chanlocs = EEG.chanlocs;

% Remove bad electrodes
EEG = pop_rejchan(EEG, 'elec',(1:size(EEG.data, 1)) ,'threshold',4,'norm','on','measure','kurt');
EEG = pop_rejchan(EEG, 'elec',(1:size(EEG.data, 1)) ,'threshold',4,'norm','on','measure','prob');        
EEG = pop_rejchan(EEG, 'elec',(1:size(EEG.data, 1)) ,'threshold',4,'norm','on','measure','spec');

% WRITE DOWN AUTOMATICALLY REMOVED ELECTRODES

eeglab redraw

% high-pass filter  
EEG = pop_eegfiltnew(EEG, [],1,1650,1,[],0);

% clean line noise
EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',(1:EEG.nbchan) ,'computepower',0,'linefreqs',[50 100] ,'normSpectrum',0,'p',0.05,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',2,'winstep',1);
                      
% plot cleaned spectrum
figure;
set(gcf, 'position', [1600 800 1000 500])
pop_spectopo(EEG, 1, [0  1915866], 'EEG' , 'percent', 35, 'freq', [6 10 22], 'freqrange',[0 100],'electrodes','on');
drawnow

% SAVE FIGURE (powerspectrum.jpg)

eeglab redraw

% SAVE FIGURE (electrodes.jpg) (Plot - channels locations - by number)
% REMOVE BAD CHANNELS (Plot - channel data (scroll) | Edit – selected data – channels – remove these)
% SAVE FIGURE (channels.png) (Plot - channel data (scroll))

%% Labels and epochs

EEG2 = EEG;

% Empty 88, 89
% Noise only 41, 42 - stim, 10 unaware, 11 weak 12 clear 
% Noise and tone, 51,52, 20, 21,22 noise, 30, 31, 32 tone

for index=1:length(EEG.event)    
    
    % ..i.e. for 1 stimulus + category task + awareness
    if strcmpi(EEG.event(index).type, '42') && (index +1 <= length(EEG.event)) && strcmpi(EEG.event(index + 1).type, '10')
        EEG.event(index).type = 'unaware_noise_noiseonly';
    elseif strcmpi(EEG.event(index).type, '42') && (index +1 <= length(EEG.event)) && strcmpi(EEG.event(index + 1).type, '11')
        EEG.event(index).type = 'weak_noise_noiseonly';
    elseif strcmpi(EEG.event(index).type, '42') && (index +1 <= length(EEG.event)) && strcmpi(EEG.event(index + 1).type, '12')
        EEG.event(index).type = 'clear_noise_noiseonly';

    elseif strcmpi(EEG.event(index).type, '52') && (index +2 <= length(EEG.event)) && strcmpi(EEG.event(index + 1).type, '20') && strcmpi(EEG.event(index + 2).type, '30')
        EEG.event(index).type = 'unaware_noise_unaware_tone';
    elseif strcmpi(EEG.event(index).type, '52') && (index +2 <= length(EEG.event)) && strcmpi(EEG.event(index + 1).type, '20') && strcmpi(EEG.event(index + 2).type, '31')
        EEG.event(index).type = 'unaware_noise_weak_tone';
    elseif strcmpi(EEG.event(index).type, '52') && (index +2 <= length(EEG.event)) && strcmpi(EEG.event(index + 1).type, '20') && strcmpi(EEG.event(index + 2).type, '32')
        EEG.event(index).type = 'unaware_noise_clear_tone';

    elseif strcmpi(EEG.event(index).type, '52') && (index +2 <= length(EEG.event)) && strcmpi(EEG.event(index + 1).type, '21') && strcmpi(EEG.event(index + 2).type, '30')
        EEG.event(index).type = 'weak_noise_unaware_tone';
    elseif strcmpi(EEG.event(index).type, '52') && (index +2 <= length(EEG.event)) && strcmpi(EEG.event(index + 1).type, '21') && strcmpi(EEG.event(index + 2).type, '31')
        EEG.event(index).type = 'weak_noise_weak_tone';
    elseif strcmpi(EEG.event(index).type, '52') && (index +2 <= length(EEG.event)) && strcmpi(EEG.event(index + 1).type, '21') && strcmpi(EEG.event(index + 2).type, '32')
        EEG.event(index).type = 'weak_noise_clear_tone';

    elseif strcmpi(EEG.event(index).type, '52') && (index +2 <= length(EEG.event)) && strcmpi(EEG.event(index + 1).type, '22') && strcmpi(EEG.event(index + 2).type, '30')
        EEG.event(index).type = 'clear_noise_unaware_tone';
    elseif strcmpi(EEG.event(index).type, '52') && (index +2 <= length(EEG.event)) && strcmpi(EEG.event(index + 1).type, '22') && strcmpi(EEG.event(index + 2).type, '31')
        EEG.event(index).type = 'clear_noise_weak_tone';
    elseif strcmpi(EEG.event(index).type, '52') && (index +2 <= length(EEG.event)) && strcmpi(EEG.event(index + 1).type, '22') && strcmpi(EEG.event(index + 2).type, '32')
        EEG.event(index).type = 'clear_noise_clear_tone';    

    elseif strcmpi(EEG.event(index).type, '89') 
        EEG.event(index).type = 'empty_center_dot';
    end    
        
end

eeglab redraw

epochs = {'unaware_noise_noiseonly', 'weak_noise_noiseonly', ...
    'clear_noise_noiseonly', 'unaware_noise_unaware_tone', ...
    'unaware_noise_weak_tone', 'unaware_noise_clear_tone', ...
    'weak_noise_unaware_tone', 'weak_noise_weak_tone', ...
    'weak_noise_clear_tone', 'clear_noise_unaware_tone', ...
    'clear_noise_weak_tone', 'clear_noise_clear_tone', ...
    'empty_center_dot'};

EEG2 = EEG; % just in case

EEG = pop_saveset(EEG, strcat('/Users/mikkohapponen/Documents/Tiedostot/Opiskelu/DI/Electroencephalography/EEG Signals/Datasets/Preprocessing/Subject', num2str(subjNumber),'BeforeBsln.set'));

EEG = pop_epoch( EEG, epochs, [-0.2, 0.8], 'epochinfo', 'yes');

EEG.subject = subjNumber;

EEG = pop_eegfiltnew(EEG, [], 27, [], false, [], 0); % 30Hz half-cutoff filter

eeglab redraw

% REMOVE BAD EPOCHS (Plot - channel data (scroll))
% SAVE FIGURE (epochs.png) (Plot by number, Plot channel data scroll)

%% ICA

% let's see how it was before ICA
pop_timtopo(EEG);

% SAVE FIGURE (beforeICA.jpg)

% ICA
[EEG] = pop_runica(EEG, 'icatype', 'runica');
 
% save before removing components
EEG = pop_saveset(EEG, strcat('/Users/mikkohapponen/Documents/Tiedostot/Opiskelu/DI/Electroencephalography/EEG Signals/Datasets/Preprocessing/Subject', num2str(subjNumber),'BeforeIca.set'));

eeglab redraw

% PICK BAD COMPONENETS BY HAND (Tools - Classify components using IClabel - Label components - Ok)
% SAVE FIGURES (IC*.png)

EEG = pop_subcomp(EEG, [8 11 16 18 23 26 27 29 41 42 45 47]); % write bad IC numbers, i.e. [1,3,5...]

eeglab redraw

% remove baseline
EEG = pop_rmbase(EEG, [-200 0]);
% interpolate removed channels
EEG = pop_interp(EEG, chanlocs, 'spherical');

% let's see how that improved
pop_timtopo(EEG);

% SAVE FIGURE (afterICA.jpg)

% pop_plottopo(EEG);

% Save ready, preprocessed data for subject
EEG = pop_saveset(EEG, strcat('/Users/mikkohapponen/Documents/Tiedostot/Opiskelu/DI/Electroencephalography/EEG Signals/Datasets/Preprocessing/Subject', num2str(subjNumber),'Ready.set'));

%% ERP analysis

% We are interested in the main contrast, NCC of awareness in MI 
mi_noise_aware_tone_aware = [];
mi_noise_aware_no_tone = []; 

m = [1,2,3,4]; % your subject numbers (datasets included in the analysis)

for i=1:length(m)
   s = m(i);
   EEG = pop_loadset(['/Users/mikkohapponen/Documents/Tiedostot/Opiskelu/DI/Electroencephalography/EEG Signals/Datasets/Preprocessing/Subject' num2str(s) 'Ready.set']);
   
   EEG1 = pop_epoch( EEG, {'weak_noise_noiseonly', 'clear_noise_noiseonly'}, [-0.2, 0.8], 'epochinfo', 'yes');
   mi_noise_aware_no_tone = cat(3, mi_noise_aware_no_tone, EEG1.data);
   
   EEG2 = pop_epoch( EEG, { 'weak_noise_weak_tone', 'weak_noise_clear_tone', ...
       'clear_noise_weak_tone', 'clear_noise_clear_tone'}, [-0.2, 0.8], 'epochinfo', 'yes');
   mi_noise_aware_tone_aware = cat(3, mi_noise_aware_tone_aware, EEG2.data);
end     

% Difference waves
diff = mean(mi_noise_aware_tone_aware,3) - mean(mi_noise_aware_no_tone,3);

% scalp topographical plot figure
times = EEG.times;
tp = [100 200 250 300 400 450 500 550 600 650 700 750 800];

tpp = round(dsearchn(times', tp'));

% this piece of code (a figure) is executed in a row. select the whole peace and then run
% start of the piece
figure('position',[680 183 1400 250]);
for t = 1:length(tp)
    x = subplot(1,length(tp),t);
    subplot(1,length(tp),t)
    topoplot(diff(:,tpp(t),:), EEG.chanlocs, 'maplimits', [-2 2], 'conv', 'off');
    title([num2str((tp(t))), ' ms']);
    set(gca,'FontSize',18);
    pos1 = get(x, 'Position'); % gives the position of current sub-plot
    new_pos1 = pos1 +[-0.005 0 0.015 0];
    set(x, 'Position',new_pos1 )
end 
% end of the piece

% SAVE FIGURE (erp1.jpg)

% Single electrode figure, let it be Fz
v = 17;  % 43 5 19 18 23 13 27 35 21 - you can change channel number
c = mean(mi_noise_aware_tone_aware, 3);
d = mean(mi_noise_aware_no_tone, 3);
% this piece of code (a figure) is executed in a row. select the whole peace and then run
% start of the piece
figure; 
p = plot(EEG.times, c(v,:), EEG.times, d(v,:), EEG.times, diff(v,:));
title(EEG.chanlocs(v).labels);
legend("MI with tone", "MI without tone", "Difference");
ylim([-2.5 2.5])
xlim([-200, 800])
xL = xlim;
yL = ylim;
line([0 0], yL);  %x-axis
line(xL, [0 0]);  %y-axis
% colormap('jet');
set(p(1),'linewidth',1);
set(p(1),'color','blue');
set(p(2),'linewidth',1);
set(p(2),'color','red');
set(p(3),'linewidth',2);
set(p(3),'color','green');
% end of the piece

% SAVE FIGURE (erp2.jpg)

%% Simple statistics

% As we are not focusing on data aggregation in matlab, don't have many 
% factors (only 1, awareness with 2 conditions), and don't have enough data, 
% we will do the t-test for the mean amplitudes in the auditory awareness 
% negativity (AAN) time window

% 1. We define AAN time window. 
% We should do it before the exp. to avoid p-hacking. Let's take the AAN tw
% from the previous experiment on AAN to be 160-300 ms
tp1 = 160;
tp2 = 300; 

% 2. We extract EEG time points of interest (remember, due to the sampling
% rate the EEG time points do not equal times they represent) 

tpp1 = round(dsearchn(times', tp1'));
tpp2 = round(dsearchn(times', tp2'));

% 3. We need to define cluster of AAN electrodes. 
% Let's use central electrodes:

cluster_labels = {'Cz', 'Pz', 'Fz', 'F1', 'F2', 'FC1', 'FC2', 'C1', 'C2', ...
                  'CPz', 'CP1', 'CP2', 'C3', 'C4', 'FC3', 'FC4', 'F3', 'F4', 'CP3', 'CP4'};

% Find indices of these electrodes in EEG.chanlocs
channelsAAN = find(ismember({EEG.chanlocs.labels}, cluster_labels));

% We average data for each paritipant (but not accross participants!) for
% these electrodes and time points

mi_noise_aware_tone_aware = [];
mi_noise_aware_no_tone = []; 

for s=1:length(m)

    t = m(s);
    % remember to put YOUR PATH to Ready files
    EEG = pop_loadset(['/Users/mikkohapponen/Documents/Tiedostot/Opiskelu/DI/Electroencephalography/EEG Signals/Datasets/Preprocessing/Subject' num2str(t) 'Ready.set']);
    
    % Re-referencing to linked mastoid reference (avg. channels TP9 and TP10) 
    % Maybe? Or maybe not? What do you think?
    EEG = pop_reref( EEG, [31,32]);

    EEG1 = pop_epoch( EEG, {'weak_noise_noiseonly', 'clear_noise_noiseonly'}, [-0.2, 0.8], 'epochinfo', 'yes');
    EEG2 = pop_epoch( EEG, { 'weak_noise_weak_tone', 'weak_noise_clear_tone', ...
               'clear_noise_weak_tone', 'clear_noise_clear_tone'}, [-0.2, 0.8], 'epochinfo', 'yes');
    
    aware = EEG2.data;
    unaware = EEG1.data; 

    % here we average to the time window and cluster avg
    aware = mean(mean(mean(aware(channelsAAN,tpp1:tpp2,:),3),2)); 
    unaware = mean(mean(mean(unaware(channelsAAN,tpp1:tpp2,:),3),2)); 

    % aware(channelsAAN, tpp1:tpp2,:) selects all trials (3rd dimension,
    % ":"), time range from tpp1 to tpp2 (2nd dimension) and channels (1rst
    % dimension)

    % mean(mean(mean....),3),2)) runs averaging 3 times: first accross 3rd
    % dimension (trials), then 2nd dimension, time points, and then the
    % remaining 1st dimension, channels. The order could be changed

    % now put those values in the array for all subs
    mi_noise_aware_no_tone = [mi_noise_aware_no_tone, unaware];
    mi_noise_aware_tone_aware = [mi_noise_aware_tone_aware, aware];

end    

[h,p,ci,stats] = ttest(mi_noise_aware_tone_aware,mi_noise_aware_no_tone)

% So, what happens?
% Was it significant? 
% Why? Share your thoughts

%% END! 

%% BONUS: FMUT

% BONUS: This will work only if you install FMUT and MUT toolbox on matlab/EEGLAB 
% binsff11.txt will be provided


for s=1:length(m) 
    t = m(s);
    EEG = bin_info2EEG(['..../Preprocessing/Subject' num2str(t) 'Ready.set'],'binsff11.txt');
    EEG1 = pop_epoch( EEG, {  'bin1', 'bin2'}, [-0.2 0.8], 'epochinfo', 'yes');
    EEG1.subject = num2str(s);
    EEG1.chanlocs = chanlocs;
    EEG1 = pop_saveset(EEG1, strcat('.../Subject', num2str(t),'Ready3.set'));
end

GND=sets2GND('gui'); % Remember to use newly saved Ready3.set this time

% Run FMUT and see
GND = FclustGND(GND, ... 
                'bins', [1, 2], ...  
                'factor_names', {'Awareness'}, ...  
                'factor_levels', [2], ...
                 'alpha', 0.05, 'chan_hood', 75, 'n_perm', 1000);




function seqbuild_single_fulltz_meas_20180829_func(pulselength, bpw, burnreps, bph, wait_time, numreads, fgen1)

    use_fgen1 = fgen1;
    params.bph = bph;
    % note that there's a lot of vesitigial code in here...
        
    % This script is set up to send in a pulse to do PL lifetime measurements. 
    % I'm using fgens 2 and 4 to send pulses and the srs fgen to trigger these.
    % 
    % I'll use the following mapping
    % fgen2 ch1: input amplitude (split and sent to the 200 MHz AOM drivers)
    % fgen2 ch2: shutter signal (will switch this to drive the EOM if necessary)
    % fgen4 ch1 and ch2 won't be doing anything important here for this
    % iteration. you could use these for frequency control in the future. The
    % variable frequency driver will be just plugged into DC 2.5V

    %% open connections to function generators (out of loop)
    % if use_fgen1
    %     try
    %         fgen1 = agilent33250a_class_new;
    %         fgen1.open;
    %     catch
    %         disp('fgen1 did something weird... trying to open again')
    %         fgen1 = agilent33250a_class_new;
    %         fgen1.open;
    %     end
    % end
    % fgen2 = agilent33522a_class_new; 
    % fgen2.open;
    
    % dgen = srsdelaygen;
    % dgen.open;

    % fgen3 = agilent33220a_class_new;
    % fgen3.open

    %fgen4 = agilent33522a_class_new;
    %fgen4.rsrcname = 'USB0::0x0957::0x2307::my50000235::0::INSTR';
    %fgen4.open

    % using an SRS function generator to trigger the experiment.
    % dgen = srsfgen;
    % dgen.rsrcname = 'GPIB1::19::INSTR';
    % dgen.open;
    
    %dgen = srsdelaygen;
    %dgen.open;

    %% create figure handles
    plotsequence = false;
    % for plotting the fgen output
    if plotsequence
        figure;
        hall = subplot(2,1,1);
        hallf = subplot(2,1,2);
    end

    %% turn off fgens
    %fgen1.output('off')
    % fgen2.ch1.output(par'off')
    %fgen4.ch1.output('off')

    %fgen4.ch2.output('off')
    pause(2)
    %% define parameters
    params.totalexptime = 5000*1e-6; % total experiment time in seconds                                                                                                                                                                                                                                    *1e-3;

    %wait parameters (s)
    params.wait = wait_time;

    params.smallestwaitstep = 1e-6;                                                   
    params.nwaitseq = round(params.wait/params.smallestwaitstep);

    %%%%
    params.timestep = 0.004*1e-6; %in seconds! % smallest it can go is 0.004
    params.srate = 1/params.timestep; %note the max sampling rate is 250 MSa/s

    % voltage levels for channels 1 and 2 on fgen2

    % ch1 is going to the input AOMs
    params.ch1amp = 2;
    params.ch1off = 0*-0.01;

    % ch2 is going to the shutter AOM
    params.ch2amp = 5; %~40MHz/V. % 2.5V~40 MHz for 80 MHz AOM
    params.ch2off = 0;%-0.135;

    % parameters for the shutter (fgen1)
    % using to frequency shift the AOMs
    params.ch5amp = 1.25;
    params.ch5off = 1.8;
    %% various time offsets
    % defined below now
    %params.shutteroffset = 1.0;%1.20; %0.18;%1.25;%0.41 %0.93

    params.ch1delay= 0; % delay of channel 1 from the SRS trigger. 
    params.ch2delay = 0; % to accomodate for time delay in AOM
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% define prep sequence. 
    %%%%% common parameters (will try to condense even further in the future)
    params.prep1time = bpw+8; % uw
    if (burnreps==0||bpw==0)
        params.nprep1seq = 1;
    else
        params.nprep1seq = burnreps;
    end

    %%%% create all sequence objects. should condense this
    prep1seq_ch1 = arbseq_class('prep1ch1',params.timestep,params.prep1time,params.nprep1seq); % 
    prep1seq_ch2 = arbseq_class('prep1ch2',params.timestep,params.prep1time,params.nprep1seq); % create sequence object
    
    %%%% ch1 - AOMs

    % prep1seq_ch1.seqtype = 'DC'; % 
    % prep1seq_ch1.heights = 0; %0.55
    prep1seq_ch1.seqtype = 'rectangular';
    prep1seq_ch1.pulseref = 'edge';
    
    %bpw = 0.5;
    
    if (burnreps==0||bpw==0)
        bpw=1;
        prep1seq_ch1.heights = 0*[1]; %0.55
    else
        prep1seq_ch1.heights = params.bph*[1]; %0.55
    end
    params.burn_pulse_width=bpw;
    prep1seq_ch1.widths = [bpw];
    prep1seq_ch1.delays = [2];

    %%%% ch2 - Shutter
    prep1seq_ch2.seqtype = 'DC';
    prep1seq_ch2.heights = 0; %0.55

    %%%% create all the sequences 
    prep1seq_ch1.createsequence;
    prep1seq_ch2.createsequence;
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Define wait sequences    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    waitseq_ch1 = arbseq_class('waitch1',params.timestep,params.smallestwaitstep/prep1seq_ch1.timeexp,params.nwaitseq);
    waitseq_ch2 = arbseq_class('waitch2',params.timestep,params.smallestwaitstep/prep1seq_ch1.timeexp,params.nwaitseq);

    %%%% ch1 - AOMS
    waitseq_ch1.seqtype = 'DC';
    waitseq_ch1.heights = 0;

    %%%% ch2 - Shutter
    waitseq_ch2.seqtype = 'DC';
    waitseq_ch2.heights =0;


    %%%% create all the sequences 
    waitseq_ch1.createsequence;
    waitseq_ch2.createsequence;
    %% define write sequence. (raman echo sequence) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    params.writetime = 20; % was 35
    params.nwriteseq = numreads;
    %%%%%%%%%%%%%%%%%%%

    writeseq_ch1= arbseq_class('writech1',params.timestep,params.writetime,params.nwriteseq); % create sequence object
    writeseq_ch2= arbseq_class('writech2',params.timestep,params.writetime,params.nwriteseq); % create pulse object
    %%%% ch1
    %writeseq.seqtype = 'DC';
    writeseq_ch1.seqtype = 'rectangular';
    writeseq_ch1.pulseref = 'edge';


    p1w = pulselength;%1.5;%2; % 150 ns is about the shortest I can go at the moment without dropping the power. 
    params.p1w = p1w;
    sw = 17; % shutter width
    writeseq_ch1.widths = [p1w];
    writeseq_ch1.delays = [1];
    writeseq_ch1.heights = [1] ; %0.55  

    %%%% ch2

    writeseq_ch2.seqtype = 'rect';
    writeseq_ch2.pulseref = 'edge'; 

    params.shutteroffset = -0.0; % 0.08 will start to show a sliver.0.05 straight thru the aom. % 0.01 works, prev -0.01
    % shutter was at -0.02 on 4/19/2018 % 0.7
    writeseq_ch2.widths  = [p1w sw]; % 0.01 0.02%[0.1 0.17] , [0.065 0.075]%100ns drop by a factor of 10 in intensity compared to 400 ns. 0.048 0.05
    writeseq_ch2.delays = [1 0]+params.shutteroffset; %0.025 previously, -.08
    writeseq_ch2.heights = [0 1]; %0.55
    
    %%%% create all sequences
    writeseq_ch1.createsequence;
    writeseq_ch2.createsequence;

    %%%% need to modify the marker settings for the write sequence
    %%%%% SHOULD MODIFY TO MATCH THE HETERODYNE SCRIPT TO KEEP THE TRIGGER
    %%%%% RIGHT BEFORE THE READ
    marker = floor((sum(writeseq_ch1.delays(1:end)) + writeseq_ch1.widths(1))/writeseq_ch1.timestep);
    %marker = floor((sum(writeseq_ch1.delays(1:end))+0.25)/writeseq_ch1.timestep);

    writeseq_ch1.repeatstring = 'repeat';
    writeseq_ch1.markerstring = 'HighAtStart';
    writeseq_ch1.markerloc = marker;

    writeseq_ch2.repeatstring = 'repeat';
    writeseq_ch2.markerstring = 'HighAtStart';
    writeseq_ch2.markerloc = marker;

    %     figure;
    %     ha1 = subplot(2,1,1);
    %     ha2 = subplot(2,1,2);
    %     writeseq_ch1.plot(ha1);
    %     writeseq_ch2.plot(ha2);
    %     ylim(ha1,[-0.1,1.1])
    %     ylim(ha2,[-0.1,1.1])

    %% "end sequence" i.e. sequence to play while waiting to trigger
    endseq_ch1 = arbseq_class('endch1',params.timestep,params.smallestwaitstep/prep1seq_ch1.timeexp);
    endseq_ch2 = arbseq_class('endch2',params.timestep,params.smallestwaitstep/prep1seq_ch1.timeexp);

    %%%% ch1
    endseq_ch1.seqtype = 'DC';
    endseq_ch1.heights = 0;

    %%%% ch2
    endseq_ch2.seqtype = 'DC';
    endseq_ch2.heights = 0;


    %%%%
    endseq_ch1.createsequence;
    endseq_ch2.createsequence;

    endseq_ch1.nrepeats = 0;
    endseq_ch2.nrepeats = 0;

    endseq_ch1.repeatstring = 'onceWaitTrig';
    endseq_ch2.repeatstring = 'onceWaitTrig';

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % stuff to send to the 80 MHz function generator.
    % note i can use a different sampling rate here b/c it's a different
    % machine!
    loseq= arbseq_class('loseq',0.04*1e-6);
    
    %loseq.seqtype = 'triangle';
    %loseq.pulseref = 'edge';
    
    loseq.seqtype = 'rectangular';
    loseq.pulseref = 'edge';


    loseq.totaltime = params.writetime;
    
    %if (blh=='burn_hi')

       loseq.widths = 0.1;
       loseq.delays = 0;%!!!! make sure this isn't zero or the fgen will just be dc at the max value
       loseq.heights = 1;
        
        %params.ch5amp = 2.8;
        %params.ch5off = params.ch5off +params.ch5amp;
    
    %elseif (blh=='burn_lo')
       % loseq.widths = 10;
       % loseq.delays = 0.1;%!!!! make sure this isn't zero or the fgen will just be dc at the max value
       % loseq.heights = 1;
    %else
    %    error('Choose burn_high or burn_low')
    %end
    %}


    % params.totalwaittime = params.shutteroffset*writeseq.timeexp + writeseq.delays(end)*writeseq.timeexp;
    % params.totalwaittime = params.shutteroffset*writeseq.timeexp + (writeseq.delays(end-1)+writeseq.delays(end))*writeseq.timeexp + 0.5*writeseq.widths(end)*writeseq.timeexp;
    params.totalwaittime =  params.nprep1seq*prep1seq_ch1.totaltime*prep1seq_ch1.timeexp+...
                            params.nwaitseq*waitseq_ch1.totaltime*waitseq_ch1.timeexp;
                             %params.nprep2seq*prep2seq_ch1.totaltime*prep2seq_ch1.timeexp;%+...



    params.fgen1delay = params.totalwaittime;
    params.fgen2delay = 0;

    loseq.nrepeats = params.nwriteseq;
    loseq.createsequence;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot the sequence
    fullseq = [repmat(prep1seq_ch1.ydata,1,params.nprep1seq),repmat(writeseq_ch1.ydata,1,params.nwriteseq), endseq_ch1.ydata];
    % fullfseq = [repmat(prep1seq_ch3.ydata,1,params.nprep1seq),repmat(writeseq_ch2.ydata,1,params.nwriteseq), endseq_ch2.ydata];
    timeseq = 0:writeseq_ch1.timestep:(length(fullseq)*writeseq_ch1.timestep-writeseq_ch1.timestep);
    % timefseq = 0:writeseq_ch1.timestep:(length(fullfseq)*writeseq_ch1.timestep-writeseq_ch1.timestep);

    if plotsequence
        plot(hall,timeseq,params.ch1off + params.ch1amp*fullseq)
        plot(hallf,timefseq,params.ch2off + 0.5*params.ch2amp*fullfseq)
        ylim(hall, [-0.1 params.ch1amp + 0.1])
        ylim(hallf, [0 5])
        xlim(hallf,[0 max(timefseq)])
        xlim(hall,[0 max(timeseq)])
    end
    tmult = 1e-6;
    fulltime = tmult*params.prep1time*params.nprep1seq + params.wait + tmult*params.writetime*params.nwriteseq;
    % if fulltime>params.totalexptime
    %     error('it is too dang long!')
    % end

    %% compile the sequence from the segments
    fullch1seq = {prep1seq_ch1,waitseq_ch1,writeseq_ch1,endseq_ch1};
    fullch2seq = {prep1seq_ch2,waitseq_ch2,writeseq_ch2,endseq_ch2};
    % fullch1seq = {writeseq_ch1,endseq_ch1};
    % fullch2seq = {writeseq_ch2,endseq_ch2};

    fgen2.ch1.clearvolatilememory
    fgen2.ch2.clearvolatilememory

    %%
    fgen2.ch1.sendandloadarbseq('ch1seq',fullch1seq)
    fgen2.ch2.sendandloadarbseq('ch2seq',fullch2seq)
    %fgen4.ch1.sendandloadarbseq('ch3seq',fullch3seq)
    %fgen4.ch2.sendandloadarbseq('ch4seq',fullch4seq)

    %% then need to send voltage and trigger settings
    % 
    fgen2.ch1.sendtrig(params.ch1delay,'EXT','POS');
    fgen2.ch2.sendtrig(params.ch2delay,'EXT','POS');
    %fgen4.ch1.sendtrig(params.ch3delay,'EXT','POS');
    %fgen4.ch2.sendtrig(params.ch4delay,'EXT','POS');

    %fgen2.send('output:sync:polarity inv')


    fgen2.ch1.sendvoltscaletoseq(fullch1seq,params.ch1amp,params.ch1off)
    fgen2.ch2.sendvoltscaletoseq(fullch2seq,params.ch2amp,params.ch2off)
    %fgen4.ch1.sendvoltscaletoseq(fullch3seq,params.ch3amp,params.ch3off)
    %fgen4.ch2.sendvoltscaletoseq(fullch4seq,params.ch4amp,params.ch4off)

    %fgen2.ch1.sendvolt;
    %fgen2.ch2.sendvolt;

    fgen2.ch1.output('on');
    fgen2.ch2.output('on');

    %fgen4.ch1.output('on');
    %fgen4.ch2.output('on');


    % dgen.send('FUNC 1') % set to square wave;
    % dgen.send('ATTL') % set to TTL levels
    % dgen.send(sprintf('FREQ %f',1/params.totalexptime));

    if use_fgen1
        try
            fgen1.sendtrig(0,'EXT','POS');
            fgen1.sendarbseq(loseq)

            fgen1.ncyc = loseq.nrepeats;
            fgen1.sendburst
            %fgen1
            %fgen1.sendvolt(0.1,0)
            fgen1.sendvolt(params.ch5amp,params.ch5off);
            pause(1)
            fgen1.output('on')        

        catch
            disp('fgen1 failed on the first attempt')
            fgen1 = agilent33250a_class_new;
            fgen1.open;
            fgen1.sendtrig(0,'EXT','POS');
            fgen1.sendarbseq(loseq)

            fgen1.ncyc = loseq.nrepeats;
            fgen1.sendburst
            %fgen1.sendvolt(0.1,0)
            fgen1.sendvolt(params.ch5amp,params.ch5off);
            pause(1)
            fgen1.output('on')        
        end
    end
    
    %% open connection to delay generator!
    try
        dgen.trig.setmode('INT');
        dgen.trig.setimp('50 ohm');
        dgen.trig.setrate(1/params.totalexptime)
        dgen.trig.setslope('RISE');
        dgen.trig.setlevel(1);

        dgen.sendtrigger;

        dgen.A.setdelay(params.fgen1delay);
        dgen.B.setdelay(params.fgen2delay);
        dgen.C.setdelay(params.fgen3delay);
        dgen.D.setdelay(params.fgen4delay);
        dgen.A.setimp('high z');
        dgen.B.setimp('high z');
        dgen.C.setimp('high z');
        dgen.D.setimp('high z');
        dgen.T0.setimp('high z');
        sendalldelays(dgen);
        sendalloutputs(dgen);
        % params.totalexptime
        disp('dgen loaded')
    catch
        disp('dgen failed on first attempt')
        dgen.trig.setmode('INT');
        dgen.trig.setimp('50 ohm');
        dgen.trig.setrate(1/params.totalexptime)
        dgen.trig.setslope('RISE');
        dgen.trig.setlevel(1);

        dgen.sendtrigger;

        dgen.A.setdelay(params.fgen1delay);
        dgen.B.setdelay(params.fgen2delay);
        dgen.C.setdelay(params.fgen3delay);
        dgen.D.setdelay(params.fgen4delay);
        dgen.A.setimp('high z');
        dgen.B.setimp('high z');
        dgen.C.setimp('high z');
        dgen.D.setimp('high z');
        dgen.T0.setimp('high z');
        sendalldelays(dgen);
        sendalloutputs(dgen);
        % params.totalexptime
        disp('dgen loaded')
    end
    
    %fgen2.disp('off')
    assignin('base','params',params);
end
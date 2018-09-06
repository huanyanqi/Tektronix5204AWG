fgen = tek5204_AWG;
fgen.open();

% Args are: Waveform name, timestep/s, total time/us
% Create a 2 rectangular pulses in a total length of 1.
thinrecpulse = arbseq_class('thinrecpulse', 0.0001, 1); 
thinrecpulse.seqtype = 'rectangular';
thinrecpulse.pulseref = 'edge';
thinrecpulse.widths = [0.1 0.1];
thinrecpulse.delays = [0.3 0.2];
thinrecpulse.heights = [1 -1];
thinrecpulse.createsequence;
% Set up some arbitrary markers that keep flipping between 0 and 1.
% Marker 2 is out of phase with marker 1.
thinrecpulse.nummarkerchannels = 2;
for i=1:length(thinrecpulse.ydata)
    if mod(floor(i/1000), 2)==0
        thinrecpulse.markerdata(1, i) = 1;
    else
        thinrecpulse.markerdata(2, i) = 1;
    end
end

% Create a fat gaussian pulse with width 1.5 in a total length of 5.
gaussianpulse = arbseq_class('gaussianpulse', 0.0001, 5);
gaussianpulse.seqtype = 'gaussian';
gaussianpulse.pulseref = 'center';
gaussianpulse.widths = [1.5];
gaussianpulse.delays = [2.5];
gaussianpulse.heights = [1];  
gaussianpulse.createsequence;

% Storing the two waveforms into the AWG waveform memory.
disp('Sending the two pulses...');
fgen.sendarbseq(thinrecpulse);
fgen.sendarbseq(gaussianpulse);
disp('Done!');

% First test to see if stuff works by just playing the waveform
% fgen.outputwav(1, thinrecpulse);
% fgen.setrunmode(1, 'CONT');
% fgen.output(1, 'ON');
% 
% fgen.outputwav(2, gaussianpulse);
% fgen.setrunmode(2, 'CONT');
% fgen.output(2, 'ON');

% Create a periodic sequence that involves doing the rectangular pulse 5 times.
% Instead of creating a proper class we just create a sequence struct.
% Format: 
% name - string indiciating the sequence's name
% seqlist - cell array of the constituent waveforms
% repistli - matrix of the number of repeats of the constituent waveforms
periodicrecseq.name = 'periodicrecseq';
periodicrecseq.seqlist = {thinrecpulse};
periodicrecseq.replist = [5];
fgen.sendandloadarbseq(periodicrecseq);

% Create a mixed sequence that involves 1 rectangular pulse, 1 gaussian pulse, and 2 rectangular pulses.
mixseq.name = 'mixseq';
mixseq.seqlist = {thinrecpulse, gaussianpulse, thinrecpulse};
mixseq.replist = [1, 1, 2];
fgen.sendandloadarbseq(mixseq);

% Set channel 1 to output the periodic rectangular and set trigger to A.
fgen.outputseq(1, periodicrecseq, 'periodicrecseq');
fgen.setseqrunmode(periodicrecseq, 'periodicrecseq', 'TRIGGERED', 'ATRigger');
fgen.output(1, 'ON');

% Set channel 2 to output the mixed sequence and set trigger to B.
fgen.outputseq(2, mixseq, 'mixseq');
fgen.setseqrunmode(mixseq, 'mixseq', 'TRIGGERED', 'BTRigger');
fgen.output(2, 'ON');

% Manually trigger the triggers.
% fgen.manualtrig('ATRigger');
% fgen.manualtrig('BTRigger');

% Turn on the AWG overall output
fgen.AWGoutput('ON')
classdef arbseq_class < handle
    % class to generate vectors/data to send to the arbitrary
    % waveform generators
    % created by J.Kindem 2016-2017
    
    properties
        seqtype = 'rectangular'; %'rectangular' or 'gaussian' pulses
        pulseref = 'edge'; % 'center' or 'edge'
        %when using rectangular pulses, whether to reference from the
        %center of the previous pulse (as done with gaussian pulses) or the
        %edge]
        
        widths = 1;
        delays = 1;
        heights = 1; %normalized height between 0 and 1.
        heightoffset = 0;
        
        autopad = false; %automatically pad the waveform if the pulse extends past the total time
        
        totaltime = 10; %total time
        timestep = 0.01;
        timeexp = 1;%1E-6; %use to get absolute time i.e. timeexp = 1e-6, time is in us.
        srate
        
        cutoff = 0.05
        normydata % in case i need to store the normalized data.
        
        symmetry = 1;
        rampstart = 0;
        rampstop = 1;
        
        name = 'arb'
        nrepeats = 0;
        repeatstring = 'repeat';
        markerstring = 'lowAtStart'; %previously, was 'highAtStartGoLow'
        markerloc = 10;
        
        seqstring

        nummarkerchannels = 0;
        % NOTE: Added for compat with Tektronix AWG. 
        % 4 x length(tdata) series of single bits for marker data.
        markerdata 
    end
    
    properties (SetAccess = 'private')
        tdata %time vector with time steps (to make things easier to plot)
        ydata %vector containing heights of arb wave
    end
    
    methods
        function obj = arbseq_class(varargin)
            switch nargin
                case 0
                    %createpulses(obj);
                case 1
                    obj.name = varargin{1};
                case 2
                    obj.name = varargin{1};
                    obj.timestep = varargin{2}/obj.timeexp;
                case 3
                    obj.name = varargin{1};
                    obj.timestep = varargin{2}/obj.timeexp;
                    obj.totaltime = varargin{3};
                    %createpulses(obj,varargin{1},varargin{2},varargin{3})
                case 4
                    obj.name = varargin{1};
                    obj.timestep = varargin{2}/obj.timeexp;
                    obj.totaltime = varargin{3};
                    obj.nrepeats = varargin{4};
            end
        end
        
        function setseqtype(obj,val)
            obj.seqtype = val;
            createpulses(obj)
        end
        
        function setpulseref(obj,val)
            obj.pulseref = val;
            createpulses(obj)
        end
        
        function setwidths(obj,val)
            obj.widths = val;
            createpulses(obj)
        end
        
        function setdelays(obj,val)
            obj.delays = val;
            createpulses(obj)
        end
        
        function setheights(obj,val)
            obj.widths = val;
            createpulses(obj)
        end
        
        function settotaltime(obj,val)
            obj.totaltime = val;
            createpulses(obj);
        end
        
        function createsequence(obj,varargin)
            switch nargin
                case 2
                    obj.widths = varargin{1};
                case 3
                    obj.widths = varargin{1};
                    obj.delays = varargin{2};
                case 4
                    obj.widths = varargin{1};
                    obj.delays = varargin{2};
                    obj.heights = varargin{3};
            end
            
            if length(obj.widths)<length(obj.delays)
                warning('oops. you need the same number of widths and delays. widths vector padded with ones');
                obj.widths = [obj.widths ones(1,length(obj.delays)-length(obj.widths))];
            elseif length(obj.widths)>length(obj.delays)
                warning('oops. you need the same number of widths and delays. delays vector padded with ones');
                obj.delays = [obj.delays ones(1,length(obj.widths)-length(obj.delays))];
            end
            
            npulses = length(obj.widths);
            
            if isempty(obj.heights)
                obj.heights  = ones(size(obj.widths));
            elseif length(obj.heights) < npulses
                obj.heights = [obj.heights ones(1,npulses- length(obj.heights))];
                warning('oops. heights vector padded with ones')
            elseif length(obj.heights) > npulses
                obj.heights = obj.heights(1:npulses);
                warning('oops. heights vector truncated')
            end
            
            
            %create pulses depending on pulse type
            switch obj.seqtype
                case {'Rect','rect','RECT','rectangular'}
                    
                    if sum(obj.delays)+sum(obj.widths)>=obj.totaltime
                        if obj.autopad
                            obj.totaltime = sum(obj.delays)+sum(obj.widths)+2;
                        else
                            warning('sequence longer than total time!')
                        end
                    end
                    
                    obj.tdata = 0:obj.timestep:(obj.totaltime-obj.timestep);
                    obj.ydata = zeros(size(obj.tdata));
                    
                    startstop = zeros(npulses,2);
                    pulseoverlap = zeros(npulses-1,1);
                    
                    if strcmp(obj.pulseref,'edge')
                        startstop(1,1) = obj.delays(1);
                        startstop(1,2) = startstop(1,1) + obj.widths(1);
                        if npulses>1
                            for pno = 2:npulses
                                startstop(pno,1) = startstop(pno-1,2)+obj.delays(pno);
                                startstop(pno,2) = startstop(pno,1)+obj.widths(pno);
                            end
                        end
                    elseif strcmp(obj.pulseref,'center')
                        for pno = 1:npulses
                            startstop(pno,1) = sum(obj.delays(1:pno))-obj.widths(pno)/2;
                            startstop(pno,2) = startstop(pno,1)+obj.widths(pno);
                        end
                    else
                        warning('oops, your pulse reference does not make sense! Check your syntax and try again.');
                    end
                    
                    % check for pulse overlap
                    if npulses>1
                        for pno = 2:npulses;
                            pulseoverlap(pno) = (startstop(pno,1)<startstop(pno-1,2));
                        end
                        overlaps = sum(pulseoverlap);
                        if overlaps>0
                            warning('OOPS! you have pulse overlap(s)!!! check your delays!')
                        end
                    end
                    
                    for pn=1:npulses
                        obj.ydata(startstop(pn,1) <= obj.tdata & obj.tdata < startstop(pn,2)) = obj.heights(pn);
                    end
                    %assignin('base','startstop',startstop)
                    
                case {'Gauss','Gaussian','gauss','gaussian'}
                    if strcmp(obj.pulseref,'edge')
                        warning('sorry, this script does not know how to deal with edge delays for gaussian pulses. using center to center instead.')
                    end
                    
                    if sum(obj.delays)>=obj.totaltime
                        if obj.autopad
                            obj.totaltime = sum(obj.delays)+2;
                        else
                            warning('sequence longer than total time!')
                        end
                    end
                    
                    obj.tdata = 0:obj.timestep:(obj.totaltime-obj.timestep);
                    obj.ydata = zeros(size(obj.tdata));
                    
                    gmult = 1/(2*sqrt(2*log(2))); %to convert width to FWHM of gaussian
                    for pn=1:npulses
                        obj.ydata = obj.ydata + obj.heights(pn)*exp(-(obj.tdata-sum(obj.delays(1:pn))).^2/(2*(gmult*obj.widths(pn))^2));
                    end
                    obj.ydata(obj.ydata<=obj.cutoff) = 0;
                    
                case {'Rectinv','rectinv','RECTINV','rectangularinverted'}
                    
                    if sum(obj.delays)+sum(obj.widths)>=obj.totaltime
                        if obj.autopad
                            obj.totaltime = sum(obj.delays)+sum(obj.widths)+2;
                        else
                            warning('sequence longer than total time!')
                        end
                    end
                    
                    obj.tdata = 0:obj.timestep:(obj.totaltime-obj.timestep);
                    obj.ydata = zeros(size(obj.tdata));
                    
                    startstop = zeros(npulses,2);
                    pulseoverlap = zeros(npulses-1,1);
                    
                    if strcmp(obj.pulseref,'edge')
                        startstop(1,1) = obj.delays(1);
                        startstop(1,2) = startstop(1,1) + obj.widths(1);
                        if npulses>1
                            for pno = 2:npulses
                                startstop(pno,1) = startstop(pno-1,2)+obj.delays(pno);
                                startstop(pno,2) = startstop(pno,1)+obj.widths(pno);
                            end
                        end
                    elseif strcmp(obj.pulseref,'center')
                        for pno = 1:npulses
                            startstop(pno,1) = sum(obj.delays(1:pno))-obj.widths(pno)/2;
                            startstop(pno,2) = startstop(pno,1)+obj.widths(pno);
                        end
                    else
                        warning('oops, your pulse reference does not make sense! Check your syntax and try again.');
                    end
                    
                    % check for pulse overlap
                    if npulses>1
                        for pno = 2:npulses;
                            pulseoverlap(pno) = (startstop(pno,1)<startstop(pno-1,2));
                        end
                        overlaps = sum(pulseoverlap);
                        if overlaps>0
                            warning('OOPS! you have pulse overlap(s)!!! check your delays!')
                        end
                    end
                    
                    for pn=1:npulses
                        obj.ydata(startstop(pn,1) <= obj.tdata & obj.tdata < startstop(pn,2)) = obj.heights(pn);
                    end
                    %assignin('base','startstop',startstop)
                    obj.ydata = 1-obj.ydata;
                    
                case {'afc','AFC','comb'} % probably want both a gaussian and rectangular comb.
                    if sum(obj.delays)>=obj.totaltime
                        if obj.autopad
                            obj.totaltime = sum(obj.delays)+2;
                        else
                            warning('sequence longer than total time!')
                        end
                    end
                    
                    obj.tdata = 0:obj.timestep:(obj.totaltime-obj.timestep);
                    obj.ydata = zeros(size(obj.tdata));
                    
                    startstop = zeros(npulses,2);
                    for pno = 1:npulses
                        startstop(pno,1) = sum(obj.delays(1:pno))-obj.widths(pno)/2;
                        startstop(pno,2) = startstop(pno,1)+obj.widths(pno);
                    end
                    
                    for pn=1:npulses
                        obj.ydata(startstop(pn,1) <= obj.tdata & obj.tdata < startstop(pn,2)) = obj.heights(pn);
                    end
                    
                    %                     gmult = 1/(2*sqrt(2*log(2))); %to convert width to FWHM of gaussian
                    %                     for pn=1:npulses
                    %                         obj.ydata = obj.ydata + obj.heights(pn)*exp(-(obj.tdata-sum(obj.delays(1:pn))).^2/(2*(gmult*obj.widths(pn))^2));
                    %                     end
                    
                    obj.ydata = 1-obj.ydata;
                    
                    
                case {'DC','dc'}
                    obj.tdata = 0:obj.timestep:(obj.totaltime-obj.timestep);
                    obj.ydata = obj.heights*ones(size(obj.tdata));
                    
                case {'sawtooth','Sawtooth','SAWTOOTH'}
                    obj.tdata = 0:obj.timestep:(obj.totaltime-obj.timestep);
                    obj.ydata = obj.heights*sawtooth((2*pi/obj.totaltime)*obj.tdata,obj.symmetry);
                    
                case {'Triangle','triangle'}
                    obj.tdata = 0:obj.timestep:(obj.totaltime-obj.timestep);
                    obj.ydata = obj.heights*sawtooth((2*pi/obj.totaltime)*(obj.tdata + 0.25*obj.totaltime),0.5) + obj.heightoffset;
                    
                case {'Triangle2','triangle2'}
                    obj.tdata = 0:obj.timestep:(obj.totaltime-obj.timestep);
                    obj.ydata = obj.heights*sawtooth((2*pi/obj.totaltime)*(obj.tdata + 0*0.25*obj.totaltime),0.5) + obj.heightoffset;
                    
                case {'Triangle3','triangle3'}
                    obj.tdata = 0:obj.timestep:(obj.totaltime-obj.timestep);
                    obj.ydata = obj.heights*(sawtooth((2*pi/obj.totaltime)*(obj.tdata + 0.5*obj.totaltime),0.5)) + obj.heightoffset;
                    
                case {'Sine','sine'}
                    obj.tdata = 0:obj.timestep:(obj.totaltime-obj.timestep);
                    obj.ydata = 0.5*obj.heights*sin((2*pi/obj.totaltime)*(obj.tdata)) + obj.heightoffset;
                    
                    
                case {'Ramp','ramp'}
                    obj.tdata = 0:obj.timestep:(obj.totaltime-obj.timestep);
                    obj.ydata = obj.rampstart + obj.tdata*(obj.rampstop-obj.rampstart)/obj.totaltime;
                    
                otherwise
                    warning('oops, your pulse type was not recognized! Check your syntax and try again.');
            end
            
            
            
            %renormalize if the data gets larger than 1. (e.g. from the
            %pulses overlapping)
            if max(abs(obj.ydata))>1
                obj.ydata = (obj.ydata-min(obj.ydata))/(max(obj.ydata)-min(obj.ydata));
            end
            
            obj.markerdata = zeros(4, length(obj.tdata));
        end
        
        function createadiabaticpulse(obj,fchirp,tfwhm,comp)
            obj.tdata = -(0.5*obj.totaltime):obj.timestep:(0.5*obj.totaltime)+0.0001;
            obj.ydata = zeros(size(obj.tdata));
            switch comp
                case 'I'
                    [obj.ydata , ~]= adiabaticpulseIQ(obj.tdata,fchirp,tfwhm);
                case 'Q'
                    [~ , obj.ydata]= adiabaticpulseIQ(obj.tdata,fchirp,tfwhm);
            end
            
        end
        
        function plot(obj,varargin)
            
            if ~isempty(varargin)
                assignin('base','varargin',varargin)
                if ishandle(varargin{1})
                    plot(varargin{1},obj.tdata,obj.ydata,varargin{2:end})
                else
                    plot(obj.tdata,obj.ydata,varargin)
                end
            else
                plot(obj.tdata,obj.ydata)
                %ylim([-1,1])
            end
        end
        
        function plotpulsestoh(h,obj,varargin)
            if ~isempty(varargin)
                plot(h,obj.tdata,obj.ydata,varargin)
            else
                plot(h,obj.tdata,obj.ydata)
            end
            
        end
        
        function y = getydata(obj)
            y = obj.ydata;
        end
        
        function t = gettdata(obj)
            t = obj.tdata;
        end
        
        function setydata(obj,y)
            obj.ydata = y;
        end
        
        function settdata(obj,t)
            obj.tdata = t;
        end
        
        function createseqstring(obj)
            %makes (part of) the string that needs to be sent to the
            %agilent waveform generator. should actually go in the class
            %def for that...
            obj.seqstring = sprintf('"%s",%g,%s,%s,%g',obj.name,obj.nrepeats,obj.repeatstring,obj.markerstring,obj.markerloc);
            %obj.seqstring = sprintf('"INT:\\%s.arb",%g,%s,%s,%g',obj.name,obj.nrepeats,obj.repeatstring,obj.markerstring,obj.markerloc);
        end
        
        function out = getseqstring(obj)
            obj.createseqstring;
            out = obj.seqstring;
        end
        
        function fixvoltageforaom(obj,voltcorrectionfn)
            %adjusts the output of ydata so that it's scaled to the AOM
            %response. i.e. so 0.5 gives 0.5 the power
            %voltcorrectionfn is the function that fixes the amplitude by
            %mapping desired output to the required voltage
            obj.createsequence; % so it doesnt' keep scaling.
            obj.ydata = voltcorrectionfn(obj.ydata)';
        end
        
        function fixvoltageforaomfreq(obj,fseq,setpoint,voltcorrectionfn,ampvsfreqfn)
            % similar to above.
            %adjusts the output of ydata so that it's scaled to the AOM
            %response. i.e. so 0.5 gives 0.5 the power
            %voltcorrectionfn is the function that fixes the amplitude by
            %mapping desired output to the required voltage
            % this one also takes
            obj.createsequence;
            if length(obj.ydata) == length(fseq.ydata)
                relerror = obj.ydata'.*setpoint./ampvsfreqfn(fseq.ydata);
                obj.ydata = voltcorrectionfn(relerror)';
            else
                error('frequency and amplitude vectors must be the same length!')
            end
        end 
    end 
end


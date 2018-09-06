classdef tek5204_AWG < handle
    properties
        id %VISA identifier
%         rsrcname = 'TCPIP0::192.168.1.202::inst0::INSTR';
        rsrcname = 'TCPIP::192.168.1.22::INSTR';
    end
    
    methods
        %% Basic Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = tek5204_AWG
            % initialiser for the object class.
        end

        function send(obj, varargin)
            % Send a message to the object. Can use printf formatting in the args.
            fprintf(obj.id, sprintf(varargin{:}), 'sync');
        end

        function out = read(obj)
            % Read a message from the object.
            out = fscanf(obj.id);
        end

        function open(obj)
            disp('Finding tek...');
            obj.id = visa('tek', obj.rsrcname);
            disp('Found tek...');

            % Set IO time out
            obj.id.Timeout = 10;
            % Calculate output buffer size
            obj.id.outputbuffersize = 8e6 + 125;
            
            disp('Opening connection...');
            fopen(obj.id);
            obj.send('*OPC?');
            obj.read();
            disp('Opened!');
            
            disp('Sending reset and clear signals...');
            % Clear registers and reset the device
            obj.send('*CLS');
            obj.send('*WAI');
            obj.send('*RST');
            obj.send('*OPC?');
            obj.read();
            disp('Finished initilizing the system');
            
            % Make it such that it will accept recommended settings
            obj.send('AWGControl:ARSettings ON')
            
            obj.errorcheck();
        end

        function close(obj)
            fclose(obj.id);
        end

        function errorcheck(obj)
            % Check for errors and print a message.
            disp('Checking for error...');
            obj.send('SYSTEM:ERROR:NEXT?');
            errorstr = obj.read();
            
            if strncmp (errorstr, '0,"No error"', 12)
                outputstr = 'No error.\n';
            else
                outputstr = ['Error reported: ', errorstr, '\n'];
            end
            fprintf (outputstr);
        end

        %% Output Control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function output(obj, channel, onoff)
            % Sets a channel output status to on or off
            % onoff is either 'ON' or 'OFF'.
            obj.send('OUTPut%d:STATe %s', channel, onoff);
        end
        
        function AWGoutput(obj, onoff)
            % Turn the overall output of the AWG on or off
            switch upper(onoff)
                case 'ON'
                    obj.send('AWGControl:RUN:IMMediate');
                case 'OFF'
                    obj.send('AWGCONTROL:STOP:IMMEDIATE');
                otherwise
                    error("Invalid AWG output option");
            end
        end

        function outputwav(obj, channel, waveform)
            % Set the built waveform into the specified channel
            obj.send('SOURce%d:CASSet:WAVeform "%s"', channel, waveform.name);
            % Set the number of marker channels based on the waveform
            obj.setmarkerchannel(channel, waveform.nummarkerchannels); 
        end

        function outputseq(obj, channel, sequence)
            % Set the built sequence into the specified channel
            % 1 means track 1 of the sequence
            obj.send('SOURce%d:CASSet:SEQuence "%s",1', channel, sequence.name);
            
            % Look at the number of channels that each constituent waveform
            % in the sequence has and then find the maximum and set that as
            % the number of marker channels.
            markerchannellist = zeros(1, length(sequence.seqlist));
            for i=1:length(sequence.seqlist)
                markerchannellist(i) = sequence.seqlist{i}.nummarkerchannels;
            end
            obj.setmarkerchannel(channel, max(markerchannellist));
        end

        %% Sequence Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function sendarbseq(obj, arbseq)
            % Takes an instance of the arbseq class and loads its waveform
            % into the AWG's memory under the instance's name.
            
            name = arbseq.name;
            fprintf('Trying to load waveform %s.\n', name)
            
            % Check if the waveform is already loaded. If so, we skip loading this.
            obj.send('WLISt:LIST?');
            wavlist = obj.read();
            % Split the returned list by commas, strip away any trailing whitespace
            % or double quotes.
            wavlist = strsplit(wavlist, ",");
            wavlist = strrep(strip(wavlist), '"', '');
            if any(strcmp(wavlist, name))
                fprintf('Waveform %s was already loaded and will be skipped.\n', name)
                return
            end
            
            % Set the waveform data to single precision
            yData = single(arbseq.ydata);
            % Sample rate
            arbseq.srate = 1 / (arbseq.timestep * arbseq.timeexp);
            
            % Make sure the y-series is a row vector, transpose otherwise
            if ~isrow(yData)
                yData = yData';
            end
            
            % This will rescale my signal if the maximum is above 1.
            % Otherwise, I think I run into problems with sending the values
            % to the machine.
            mx = max(abs(yData));
            if mx > 1
                warning('Data was rescaled b/c the absolute value of your waveform went above 1!')
                yData = yData / mx;
            end
            
            % Convert datapoints to binary before sending
            binblockBytes = typecast(yData, 'uint8');  
            % Number of datapoints
            numPoints = length(yData);
            % Number of bytes. Multiply by 4 because single precision float takes up 4 bytes.
            numBytes = numPoints * 4; 

            % Create a new waveform with the desired name and number of points
            obj.send('WLISt:WAVeform:NEW "%s",%d', name, numPoints)
            % Zero refers to the start index of the data, used for chunking large waveforms
            header = sprintf('WLISt:WAVeform:DATA "%s",0,%d,#%d%d', name, numPoints, length(num2str(numBytes)), numBytes);
            % Write in this weird way to send float data in binary
            fwrite(obj.id, [header, binblockBytes], 'uint8'); 
            % Blocking command, make sure the previous command is complete before proceeding.
            obj.send('*WAI');
           
            % Marker array of single bits. For the markers, each marker data occupies 
            % one bit. Four most significant bits of each byte are used for markers. 
            % Bit 7 for marker 1, bit 6 for marker 2, etc. We first convert the 
            % 0 and 1's to uint8 then shift it to the respective bit.
            markerDataArr = uint8(arbseq.markerdata);
            markerData = bitsll(markerDataArr(1,:), 7) + bitsll(markerDataArr(2,:), 6) + bitsll(markerDataArr(3,:), 5) + bitsll(markerDataArr(4,:), 4);
            % Write the marker array to this waveform, 0 is the start index
            marker_header = sprintf('WLISt:WAVeform:MARKer:DATA "%s",0,%d,#%d%d', name, numPoints, length(num2str(numPoints)), numPoints);
            fwrite(obj.id, [marker_header, markerData], 'uint8');
            obj.send('*WAI');

            % Set sampling rate
            obj.send('WLISt:WAVeform:SRATe "%s",%d', name, arbseq.srate)
            % Print waveform has been downloaded
            fprintf('Arb waveform "%s" downloaded.\n', name) 

            % Check for error
            obj.errorcheck();          
        end

        function buildarbseq(obj, sequence)
            % Takes the sequence name and a list of the actual sequences 
            % and builds a sequence in the AWG memory that combines all of that.

            nSegments = length(sequence.seqlist);

            % Create a new sequence with the specified name and number of steps
            % The last 1 means we're creating a sequence with only 1 track
            obj.send('SLISt:SEQuence:NEW "%s",%d,1', sequence.name, nSegments)

            for segment = 1:nSegments
                % Assign the segment-th waveform to the segment-th step
                % The 1 in the TASSet means we're assigning to track 1 of the sequence
                obj.send('SLISt:SEQuence:STEP%d:TASSet1:WAVeform "%s","%s"', segment, sequence.name, sequence.seqlist{segment}.name)
                % Assign the number of repeats to the segment-th stepFGEN
                obj.send('SLISt:SEQuence:STEP%d:RCOunt "%s",%d', segment, sequence.name, sequence.replist(segment))
            end
        end

        function sendandloadarbseq(obj, sequence)
            % Takes a list of waveforms and stores each one in memory, then 
            % builds a new sequence that uses these waveforms by repeating each
            % one a certain number of times.

            for currentseq = 1:length(sequence.seqlist)
                obj.sendarbseq(sequence.seqlist{currentseq});
            end
            
            % Build a string telling the AWG how to combine and repeat those 
            % individual sequences.
            obj.buildarbseq(sequence);
            
            % Look at the sampling rate for each waveform for the sequence and take the maximum rate
            % Give a warning if waveforms don't have the same rate.
            srateseq = cell2mat(cellfun(@(c) c.srate, sequence.seqlist, 'UniformOutput', false));
            if ~all(srateseq == srateseq(1))
                warning('Some waveforms in the sequence have different rates! Taking the fastest one.')
            end
            srateseq = max(srateseq);

            obj.send('SLISt:SEQuence:SRATe "%s",%d', sequence.name, srateseq)

            % Check for error
            obj.errorcheck();
        end

        function clearvolatilememory(obj, channel)
            % Clear the waveform/sequence currently assigned to that channel.
            obj.send('SOURce%d:CASSet:CLEar', channel)
        end

        function deleteallwaveforms(obj)
            % Delete all saved waveforms from memory.
            obj.send('WLISt:WAVeform:DELete ALL')
        end

        %% Sequence/Waveform Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function setrunmode(obj, channel, rmode)
            % Set the running mode of the channel to the given mode.
            % Only works for outputting waveform and not sequence.
            % rmode = CONTinuous, TRIGgered, TCONtinuous, GATed
            obj.send('SOURce%d:RMODe %s', channel, rmode);
        end
        
        function setseqrunmode(obj, sequence, rmode, varargin)
            % Sets the running and triggering mode for a particular sequence 
            % by assigning  triggers/jumps to the first and last steps.

            switch upper(rmode)
                case {'TRIG', 'TRIGGERED'}
                    trigsource = varargin{1};
                    % Set the 1st step to be triggered by the given trigger source
                    % trigsource: 'ATRigger', 'BTRigger', 'ITRigger', 'OFF'
                    obj.send('SLISt:SEQuence:STEP1:WINPut "%s",%s', sequence.name, trigsource);
                    % Set the last step to jump back to the 1st step and wait for trigger
                    obj.send('SLISt:SEQuence:STEP%d:GOTO "%s",1', length(sequence.seqlist), sequence.name);
                case {'CONT', 'CONTINUOUS'}
                    % Set the first step to not wait for anything\
                    obj.send('SLISt:SEQuence:STEP1:WINPut "%s",OFF', sequence.name);
                    % Set the last step to jump back to the 1st step and immediately repeat
                    obj.send('SLISt:SEQuence:STEP%d:GOTO "%s",1', length(sequence.seqlist), sequence.name);
                otherwise
                    error('Invalid running mode!')
            end
        end

        function manualtrig(obj, trigsource)
            % Manually send a trigger signal for testing purposes.
            % trigchannel should be either 'ATRigger' or 'BTRigger'.
            obj.send('TRIGger:IMMediate %s', trigsource)
        end

        function internaltriginterval(obj, interval)
            % Sets the internal trigger interval. Range is 1 us to 10 s.
            obj.send('TRIGger:INTerval %g', interval)
        end
        
        function setmarkerchannel(obj, channel, nummarkerchannels)
            % Set the number of bits of resolution to be used for the particular
            % channel. 16 = no marker, 15 = 1 marker, ..., 12 = 4 markers.
            obj.send('SOURce%d:DAC:RESolution %d', channel, 16-nummarkerchannels)
        end
        
        % function sendvoltscaletoseq(obj, channel, seqlist, varargin)

        %     switch nargin
        %         case 3
        %             obj.amp = varargin{1};
        %         case 4
        %             obj.amp = varargin{1};
        %             obj.off = varargin{2};
        %     end

        %     tempseq = [];
        %     for currentseq = 1:length(seqlist)
        %         tempseq = [tempseq, seqlist{currentseq}.ydata];
        %     end
        %     obj.amp = max(abs(tempseq)) * obj.amp;
            
        %     obj.send('SOURce%d:VOLT:UNIT VPP', channel);
        %     obj.send('SOURce%d:VOLT %8.5G', channel, obj.amp);
        %     obj.send('SOURce%d:VOLT:OFFS %8.5G', channel, obj.off);

        %     [SOURce[n]:]VOLTage[:LEVel][:IMMediate]:LOW
        %     [SOURce[n]:]VOLTage[:LEVel][:IMMediate]:HIGH
        %     [SOURce[n]:]VOLTage[:LEVel][:IMMediate][:AMPLitude]
        %     [SOURce[n]:]VOLTage[:LEVel][:IMMediate]:OFFSet % For DC
        %     [SOURce[n]:]VOLTage[:LEVel][:IMMediate]:BIAS % For AC
        % end

    end
end
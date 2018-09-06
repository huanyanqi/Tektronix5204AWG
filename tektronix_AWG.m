classdef tektronix_AWG < handle
   %% class def for the agilent 33522a. _new updates to use pulse_classdef
   % and actually has the capability to load pulses.
   % written by J Kindem 2014-2017
    
    %todo:
    % move commands to the channel class. 
    % commands here should send both channels...
    properties
        id %VISA identifier
        rsrcname = 'TCPIP0::192.168.1.202::inst0::INSTR';
    end
    
    methods
         function obj = tektronix_AWG
         end
             
        function open(obj)
            obj.id = instrfind('Type', 'visa-tcpip', 'rsrcname', obj.rsrcname, 'Tag', '');
            % Create the GPIB object if it does not exist; otherwise use the object that was found.

            if isempty(obj.id)
              obj.id = visa('AGILENT', obj.rsrcname);
            else
              fclose(obj.id);
              obj.id = obj.id(1);
            end

            obj.id.Timeout = 15; %set IO time out
            %calculate output buffer size
            buffer = 8e6;
            obj.id.outputbuffersize = buffer + 125;
            fopen(obj.id);
        end

        function sawtooth(obj, channel, freq, power, pol)
            % FREQUENCY IS IN MHz!
            obj.send(sprintf('FGEN:CHAN%1.0f:TYPE %s', channel, 'TRI'));

            switch pol
                case 'up'
                    obj.send(sprintf('FGEN:CHAN%1.0f:SYMM %.5f', channel, 100));
                case 'down'
                    obj.send(sprintf('FGEN:CHAN%1.0f:SYMM %.5f', channel, 0));
            end
            obj.send(sprintf('FGEN:CHAN%1.0f:FREQ %.6fe6', channel, freq));
            obj.send(sprintf('FGEN:CHAN%1.0f:AMPL:POW %.5f', channel, power));
        end
        
        function output(obj,channel,onoff)
            obj.send(sprintf('OUTP%g:STATE %s', channel, onoff))
        end

        function close(obj)
            fclose(obj.id);
        end

        function send(obj, string)
            fprintf(obj.id, string);
        end
        
    end
end


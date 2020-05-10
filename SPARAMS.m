classdef SPARAMS < handle
    % Author: Nathan Chordas-Ewell
    % Written: 2020
    %Works on up to 4 port measurements in RI, MA, dBA formats (Som
    %Usage:
    % s = SPARAMS('filename')
    % s.setZ0(Z0)
    % If you do not want to associate a file
    %  s = SPARAMS()
    %  s.setNumPorts(N)
    %  s.f = freq; s.S11 = s11; ....
    %Copy an object to a new object (without referencing a pointer)
    % snew = s.copyobj();
    %Get some value:
    % s11 = s.S11 in R + jI format
    %Plot any trace in dB:
    % s.plotdB('S11')
    % or for multiple, put arguments in cell: s.plotdB({'S11', 'S21'})
    %Plot all:
    % s.plotAlldB
    %To get the S matrix at a specific frequency:
    % s.toMat(freq)
    %To renormalize to a different system impedance
    % s.renorm(Z0new)
    % This will also automatically reset Z0 to the new value
    %To cascade N times
    % s.cascade(N)
    %Plot dispersion diagram when s is a unit cell
    % s.plotDispersion
    %Export a csv file of the S parameter data
    % s.writeCSV('filename.csv')
    %Deembed S parameters to the DUT. Returns an SPARAMS object of the DUT
    % P1 <--err1--DUT--err2--> P2
    % DUT = SPARAMS.deembed(err1, measured, err2)

    properties
        f;
        S11;S12;S13;S14;S21;S22;S23;S24;S31;S32;S33;S34;S41;S42;S43;S44;
        Z0;
        numPorts;
    end
    properties(Access = private) 
        sparams;
        rawData; data; form; filename; formStr;
        freqScale;
    end

    methods
        function obj = SPARAMS(filename)
            switch nargin
                case 1
                    obj.setFile(filename);
                case 0
                    warning('No data set. Once data is assigned you will need to call setNumPorts() to set the number of ports.');
            end
            obj.setFreqUnits('Hz');
        end
        
        function setFile(obj, fname)
            %If you wish to set a filename after creating the object. Also
            %called when initializing with a sNp file.
            obj.filename = fname;
            obj.txt2data();
            obj.getFormat();
            obj.parseData();
        end
        
        function s2 = copyobj(obj)
            %Copy an object without creating a pointer
            % s2 = s2.copyobj();
            %Note that this does not copy items such as the original
            %filename, raw data, and original s2p file format
            s2 = eval(class(obj));
            for p = properties(obj).'
                try
                    s2.(p{1}) = obj.(p{1});
                catch
                    warning('Failed to copy property: %s', p);
                end
            end
        end
        
        function setNumPorts(obj, numPorts)
            %To set the number of ports if object is created without a file
            %s.setNumPorts(N);
            %Where N = 1, 2, 3, or 4
            obj.numPorts = numPorts;
        end
        
        function setZ0(obj, sysZ0)
            %Set system Z0 if file does not have the information or object
            %is not created with a file
            %s.setZ0(z0)
            %Where Z0 usually = 50, 75
            obj.Z0 = sysZ0;
        end
        
        function setFreqUnits(obj, scale)
            %Sets the frequency scale (Hz, kHz, MHz, GHz) for plotting.
            %Does not change any S-parameter or frequency data
            %s.setFreqUnits('GHz')
           if strcmpi(scale, 'GHz')
               obj.freqScale = {10^9, 'GHz'};
           elseif strcmpi(scale, 'MHz')
               obj.freqScale = {10^6, 'MHz'};
           elseif strcmpi(scale, 'kHz')
               obj.freqScale = {10^3, 'kHz'};
           elseif strcmpi(scale, 'Hz')
               obj.freqScale = {1, 'Hz'};
           else
               warning('Frequency scale not set. Defaulting to Hz');
           end
        end
        
        function txt2data(obj)
			[~, ~, fExt] = fileparts(obj.filename);
			obj.numPorts = str2num(fExt(3));
            txt = fileread(obj.filename);
            key = '#';
            
            txt = regexprep(txt, '![\s\S]*?!|([^:]|^)!.*$', '', 'lineanchors', 'dotexceptnewline');
			txt = regexprep(txt, '\t+', ' '); %Remove all tabs and replace with s single space (dealt with in next line)
            txt = strtrim(txt);	%Remove trailing whitespace

            hashIndex = strfind(txt,key);
            newLineIndex = regexp(txt(hashIndex:end), '[\r]');
            newLineIndex = newLineIndex(1);
            obj.form = txt(hashIndex:hashIndex+newLineIndex);
            
            dataTxt = txt(hashIndex+newLineIndex:end);
			%Remove comments (!)
			dataTxt = regexprep(dataTxt, '![\s\S]*?!|([^:]|^)!.*$', '', 'lineanchors', 'dotexceptnewline');
			dataTxt = regexprep(dataTxt, '\t+', ' '); %Remove all tabs and replace with s single space (dealt with in next line
            dataTxt = strtrim(dataTxt);	%Remove trailing whitespace
            for j = 1:100 %Not a very efficient way of doing it, but I'm getting tired of this shit
				len = 101 - j;
				for k = 1:len
					regStr(k) = ' ';
				end
        		dataTxt = regexprep(dataTxt, regStr, ',');
				clear regStr
            end
            obj.rawData = cell2mat(textscan(dataTxt, '', 'delimiter', ',', 'headerlines', 0, 'emptyvalue', nan, 'collectoutput', 1));
			obj.formatMatrix();
            obj.f = obj.data(:, 1);
            if contains(obj.form, 'MHz', 'IgnoreCase', true)
                obj.f = obj.f.*1e6;
            elseif contains(obj.form, 'GHz', 'IgnoreCase', true)
                obj.f = obj.f.*1e9;
            elseif contains(obj.form, 'kHz', 'IgnoreCase', true)
                obj.f = obj.f.*1e3;
            end
            obj.sparams = obj.data(:, 2:end);
        end

		function formatMatrix(obj)
			[rows, ~] = size(obj.rawData);
			if(obj.numPorts == 1 || obj.numPorts == 2)
					obj.data = obj.rawData;
			elseif(obj.numPorts == 3)
				%Do the stuff for s3p
				rr = 1;
				for r = 1:3:rows-2
					obj.data(rr, :) = [obj.rawData(r, 1:end) obj.rawData(r+1, 1:end-1) obj.rawData(r+2, 1:end-1)];
					rr = rr + 1;
				end
			elseif(obj.numPorts == 4)
				rr = 1;
				for r = 1:4:rows-3
					obj.data(rr, :) = [obj.rawData(r, 1:end-1) obj.rawData(r+1, 2:end-1) obj.rawData(r+2, 2:end-1) obj.rawData(r+3, 2:end-1)];
					rr = rr + 1;
				end
			end

		end

        function getFormat(obj)
            obj.formStr = string(obj.form);
            if contains(obj.formStr, 'f', 'IgnoreCase', true) == 0 && contains(obj.formStr, 'Hz', 'IgnoreCase', true) == 0
                error('Looks like this is a time domain file, or not formatted properly.')
            end
            if(contains(obj.formStr, 'R', 'IgnoreCase', true))
                RIndex = strfind(obj.formStr, 'R');
                tempZ0str = char(obj.formStr);
                tempZ0str = strtrim(tempZ0str);
                Z0temp = str2double(tempZ0str(RIndex+2:end));
                obj.setZ0(Z0temp);
            else
                warning('System Z0 not set. You may do so manually with obj.setZ0()')
            end
            if contains(obj.formStr, 'dB', 'IgnoreCase', true)
                obj.form = 'dB';
            elseif contains(obj.formStr, 'MA', 'IgnoreCase', true)
                obj.form = 'MA';
            elseif contains(obj.formStr, 'RI', 'IgnoreCase', true)
                obj.form = 'RI';
            else
                error('Format unknown.');
            end
        end

        function S = formatCorrectly(obj, c1, c2)
            %Make S parameter in RI no matter orignal format
            switch obj.form
                case "MA"
                    S = c1.*exp(j*deg2rad(c2));
                case "RI"
                    S = c1 + j*c2;
                case "dB"
                    S = 10.^(c1/20).*exp(j*deg2rad(c2));
            end
        end

        function parseData(obj)
            switch obj.numPorts
                case 1
                    obj.S11 = obj.formatCorrectly(obj.sparams(:, 1), obj.sparams(:, 1+1));
                case 2
                    obj.S11 = obj.formatCorrectly(obj.sparams(:, 1), obj.sparams(:, 2));
                    obj.S21 = obj.formatCorrectly(obj.sparams(:, 3), obj.sparams(:, 4));
                    obj.S12 = obj.formatCorrectly(obj.sparams(:, 5), obj.sparams(:, 6));
                    obj.S22 = obj.formatCorrectly(obj.sparams(:, 7), obj.sparams(:, 8));
                case 3
                    obj.S11 = obj.formatCorrectly(obj.sparams(:, 1), obj.sparams(:, 2));
                    obj.S12 = obj.formatCorrectly(obj.sparams(:, 3), obj.sparams(:, 4));
                    obj.S13 = obj.formatCorrectly(obj.sparams(:, 5), obj.sparams(:, 6));
                    obj.S21 = obj.formatCorrectly(obj.sparams(:, 7), obj.sparams(:, 8));
                    obj.S22 = obj.formatCorrectly(obj.sparams(:, 9), obj.sparams(:, 10));
                    obj.S23 = obj.formatCorrectly(obj.sparams(:, 11), obj.sparams(:, 12));
                    obj.S31 = obj.formatCorrectly(obj.sparams(:, 13), obj.sparams(:, 14));
                    obj.S32 = obj.formatCorrectly(obj.sparams(:, 15), obj.sparams(:, 16));
                    obj.S33 = obj.formatCorrectly(obj.sparams(:, 17), obj.sparams(:, 18));
                case 4
                    obj.S11 = obj.formatCorrectly(obj.sparams(:, 1), obj.sparams(:, 2));
                    obj.S12 = obj.formatCorrectly(obj.sparams(:, 3), obj.sparams(:, 4));
                    obj.S13 = obj.formatCorrectly(obj.sparams(:, 5), obj.sparams(:, 6));
                    obj.S14 = obj.formatCorrectly(obj.sparams(:, 7), obj.sparams(:, 8));
                    obj.S21 = obj.formatCorrectly(obj.sparams(:, 9), obj.sparams(:, 10));
                    obj.S22 = obj.formatCorrectly(obj.sparams(:, 11), obj.sparams(:, 12));
                    obj.S23 = obj.formatCorrectly(obj.sparams(:, 13), obj.sparams(:, 14));
                    obj.S24 = obj.formatCorrectly(obj.sparams(:, 15), obj.sparams(:, 16));
                    obj.S31 = obj.formatCorrectly(obj.sparams(:, 17), obj.sparams(:, 18));
                    obj.S32 = obj.formatCorrectly(obj.sparams(:, 19), obj.sparams(:, 20));
                    obj.S33 = obj.formatCorrectly(obj.sparams(:, 21), obj.sparams(:, 22));
                    obj.S34 = obj.formatCorrectly(obj.sparams(:, 23), obj.sparams(:, 24));
                    obj.S41 = obj.formatCorrectly(obj.sparams(:, 25), obj.sparams(:, 26));
                    obj.S42 = obj.formatCorrectly(obj.sparams(:, 27), obj.sparams(:, 28));
                    obj.S43 = obj.formatCorrectly(obj.sparams(:, 29), obj.sparams(:, 30));
                    obj.S44 = obj.formatCorrectly(obj.sparams(:, 31), obj.sparams(:, 32));
            end
        end

        function [M, fAct] = toMat(obj, freq)
            %Gives S matrix at a specified frequency.
            %If frequency is not found, gives data at closest frequency
            %[S, fActual] = s.toMat(freq);
            %Where S is a numPort x numPort matrix at the specified
            %frequency (or the closest)
            [~, ind] = min(abs(obj.f-freq));
			fAct = obj.f(ind);
            if obj.f(ind) ~= freq
                warning('Not a frequency in list, using %d', obj.f(ind));
            end
            switch obj.numPorts
                case 1
                    M = obj.S11(ind);
                case 2
                    M = [obj.S11(ind), obj.S12(ind);
                        obj.S21(ind), obj.S22(ind)];
                case 3
                    M = [obj.S11(ind), obj.S12(ind), obj.S13(ind);
                        obj.S21(ind), obj.S22(ind), obj.S23(ind);
                        obj.S31(ind), obj.S32(ind), obj.S33(ind)];
                case 4
                    M = [obj.S11(ind) obj.S12(ind) obj.S13(ind) obj.S14(ind);
                        obj.S21(ind) obj.S22(ind) obj.S23(ind) obj.S24(ind);
                        obj.S31(ind) obj.S32(ind) obj.S33(ind) obj.S34(ind);
                        obj.S41(ind) obj.S42(ind) obj.S43(ind) obj.S44(ind)];
            end
        end
        
        function Sxx = todB(obj, s)
           Sxx = 20*log10(abs(eval(['obj.' s])));
        end

        function plotAlldB(obj)
            %Plot all S-parameters in dB format
            %s.plotAlldB();
           hold on
           switch obj.numPorts
               case 1
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S11)));
                   leg = {'S11'};
               case 2
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S11)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S12)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S21)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S22)));
                   leg = {'S11', 'S12', 'S21', 'S22'};
               case 3
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S11)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S12)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S13)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S21)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S22)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S23)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S31)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S32)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S33)));
                   leg = {'S11', 'S12', 'S13', 'S21', 'S22', 'S23', 'S31', 'S32', 'S33'};
               case 4
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S11)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S12)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S13)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S14)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S21)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S22)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S23)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S24)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S31)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S32)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S33)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S34)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S41)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S42)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S43)));
                   plot(obj.f/obj.freqScale{1}, 20*log10(abs(obj.S44)));
                   leg = {'S11','S12','S13','S14','S21','S22','S23','S24','S31','S32','S33','S34','S41','S42','S43','S44'};
           end
           xlabel(['f (' obj.freqScale{2} ')']);
           ylabel('S-parameters (dB)');
           legend(leg);
        end

        function plotdB(obj, splt)
            %Plots any speciified S-parameter(s) in dB format
            %Say you wish to plot only S11 and S21
            %s.plotdB({'S11' 'S21'})
            hold on
            switch class(splt)
                case 'cell'
                    leg = cell(size(length(splt), 1));
                    for i = 1:length(splt)
                        stringToPlot = ['obj.' splt{i}];
                        StoPlot = eval(stringToPlot);
						if isempty(StoPlot)
							error('%s not found', stringToPlot)
						end
                        plot(obj.f/obj.freqScale{1}, 20*log10(abs(StoPlot)));
                        leg{i} = splt{i};
                    end
                    legend(leg);
                case 'char'
                    stringToPlot = ['obj.' splt];
                    StoPlot = eval(stringToPlot);
					if isempty(StoPlot)
						error('%s not found', stringToPlot)
					end
                    plot(obj.f/obj.freqScale{1}, 20*log10(abs(StoPlot)));
            end
            xlabel(['f (' obj.freqScale{2} ')']);
            ylabel('dB');

        end

        function plotPolar(obj, splt)
            %Plots any speciified S-parameter(s) in polar format
            %Say you wish to plot only S11 and S21
            %s.plotPolar({'S11' 'S21'})
            switch class(splt)
                case 'cell'
                    leg = cell(size(length(splt), 1));
                    for i = 1:length(splt)
                        stringToPlot = ['obj.' splt{i}];
                        StoPlot = eval(stringToPlot);
                        polarplot(angle(StoPlot), abs(StoPlot));
						hold on
                        leg{i} = splt{i};
                    end
                    legend(leg);
                case 'char'
                    stringToPlot = ['obj.' splt];
                    StoPlot = eval(stringToPlot);
                    polarplot(angle(StoPlot), abs(StoPlot));
            end
        end
        
        function plotSmith(obj, splt)
            %Plots any speciified S-parameter(s) in polar format
            %Say you wish to plot only S11 and S21
            %s.plotPolar({'S11' 'S21'})
            
            % Draw Smith chart axes
            % Draw outer circle
            t = linspace(0, 2*pi, 100);
            x = cos(t);
            y = sin(t);
            plot(x, y, 'Color','black'); axis equal; 
            % Place title and remove ticks from axes
            title(' Smith Chart ')
            set(gca,'xticklabel',{[]});
            set(gca,'yticklabel',{[]});
            hold on 

            % Draw circles along horizontal axis
            k = [.25 .5 .75 ];
            for i = 1 : length(k)
                x(i,:) = k(i) + (1 - k(i)) * cos(t);
                y(i,:) = (1 - k(i)) * sin(t);
                plot(x(i,:), y(i,:), 'k')
            end 

            line ([-1 1],[0 0], 'color', 'black')

            % Draw partial circles along vertical axis
            kt = [0 2.5 pi 3.79 4.22];
            k = [0 .5  1 2 4 ];
            for i = 1 : length(kt)
                t = linspace(kt(i), 1.5*pi, 50);
                a(i,:) = 1 + k(i) * cos(t);
                b(i,:) = k(i) + k(i) * sin(t);
                plot(a(i,:), b(i,:),'k', a(i,:), -b(i,:),'k' )
            end
            hold on
            switch class(splt)
                case 'cell'
                    leg = cell(size(length(splt), 1));
                    for i = 1:length(splt)
                        stringToPlot = ['obj.' splt{i}];
                        StoPlot = eval(stringToPlot);
                        [x, y] = pol2cart(angle(StoPlot), abs(StoPlot));
                        plot(x, y);
						hold on
                        leg{i} = splt{i};
                    end
                    legend(leg);
                case 'char'
                    stringToPlot = ['obj.' splt];
                    StoPlot = eval(stringToPlot);
                    [x, y] = pol2cart(angle(StoPlot), abs(StoPlot));
                    plot(x, y);
            end
            xlim([-1.05 1.05]); ylim([-1.05 1.05])
        end

		function [Zin1, Zin2] = toInputImpedance(obj)
            %Converts to port 1 and port 2 input impedance (assuming other
            %port is terminated in Z0)
            %[Zin1, Zin2] = s.toInputImpedance();
			Zin1 = obj.Z0.*(1 + obj.S11)./(1 - obj.S11);
			Zin2 = obj.Z0.*(1 + obj.S22)./(1 - obj.S22);
		end

		function [A, B, C, D] = toABCD(obj, freq)
            %Calculates the ABCD parameters from S-parameters
            %[A, B, C, D] = s.toABCD(); for all frequencies (s.f)
            %[A, B, C, D] = s.toABCD(freq); for a specified frequency
			if obj.numPorts ~= 2
				warning('Only for 2 port measurements');
			end
			switch nargin
				case 1
					fS11 = obj.S11; fS12 = obj.S12; fS21 = obj.S21; fS22 = obj.S22;
				case 2
					S = obj.toMat(freq);
					fS11 = S(1, 1); fS12 = S(1, 2); fS21 = S(2, 1); fS22 = S(2, 2);
			end
			A = ((1 + fS11).*(1 - fS22) + fS12.*fS21)./(2.*fS21);
			B = obj.Z0*((1 + fS11).*(1 + fS22) - fS12.*fS21)./(2.*fS21);
			C = (1/obj.Z0)*((1 - fS11).*(1 - fS22) - fS12.*fS21)./(2*fS21);
			D = ((1 - fS11).*(1 + fS22) + fS12.*fS21)./(2*fS21);
        end
        
        function [T11, T12, T21, T22] = toTParams(obj, freq)
            %Calculates T parametres from S-parameters
            %[T11, T12, T21, T22] = s.toTParams(); for all frequencies
            %(s.f)
            %[T11, T12, T21, T22] = s.toTParams(freq); for a specified frequency
            switch nargin
				case 1
					fS11 = obj.S11; fS12 = obj.S12; fS21 = obj.S21; fS22 = obj.S22;
				case 2
					S = obj.toMat(freq);
					fS11 = S(1, 1); fS12 = S(1, 2); fS21 = S(2, 1); fS22 = S(2, 2);
            end
            T11 = -(fS11.*fS22 - fS12.*fS21)./fS21;
            T12 = fS11./fS21;
            T21 = -fS22./fS21;
            T22 = 1./fS21;
        end

		function renorm(obj, Z0new)
            %Renormalize to a new characteristic impedance
            %s.renorm(newZ0);
            %http://qucs.sourceforge.net/tech/node98.html
            %NOT FINISHED FOR 3 OR 4 PORTS
            Z0old = obj.Z0;
            for i = 1:length(obj.f)
				S = obj.toMat(obj.f(i));
                R = ((Z0new - Z0old)./(Z0new + Z0old)).*eye(obj.numPorts);
                A = (sqrt(Z0new/Z0old)./(Z0new + Z0old)).*eye(obj.numPorts);
                Sn = A^(-1) * (S - R) * (eye(obj.numPorts) - R * S)^(-1) * A;
% 				Sn = s2s(S, Z0old, Z0new);
                switch obj.numPorts
					case 1
						S11n(i) = Sn(1, 1);
					case 2
						S11n(i) = Sn(1, 1); S12n(i) = Sn(1, 2); S21n(i) = Sn(2, 1); S22n(i) = Sn(2, 2);
                end
                obj.Z0 = Z0new;
            end
            switch obj.numPorts
				case 1
					obj.S11 = S11n;
				case 2
					obj.S11 = S11n; obj.S12 = S12n; obj.S21 = S21n; obj.S22 = S22n; 
            end
            obj.setZ0(Z0new);
		end
		
		function [Z11, Z12, Z21, Z22] = toZparams(obj, freq)
            %Calculate Z parameters from S-parameters
            %[Z11, Z12, Z21, Z22] = toZparams(); for all frequencies (s.f)
            %[Z11, Z12, Z21, Z22] = toZparams(freq); for a specified frequency
			if obj.numPorts ~= 2
				warning('Only for 2 port measurements');
			end
			switch nargin
				case 2
					fS11 = obj.S11; fS12 = obj.S12; fS21 = obj.S21; fS22 = obj.S22;
				case 3
					S = obj.toMat(freq);
					fS11 = S(1, 1); fS12 = S(1, 2); fS21 = S(2, 1); fS22 = S(2, 2);
			end
			Z11 = obj.Z0*((1 + fS11).*(1 - fS22) + fS12.*fS21)./((1 - fS11).*(1 - fS22) - fS12.*fS21);
			Z12 = obj.Z0*(2.*fS12)./((1 - fS11).*(1 - fS22) - fS12.*fS21);
			Z21 = obj.Z0*(2*fS21)./((1 - fS11).*(1 - fS22) - fS12.*fS21);
			Z22 = obj.Z0*((1 - fS11).*(1 + fS22) + fS12.*fS21)./((1 - fS11).*(1 - fS22) - fS12.*fS21);
        end

        function [Y11, Y12, Y21, Y22] = toYparams(obj)
            %Calculate Y parameters from S-parameters
            %[Y11, Y12, Y21, Y22] = toYparams(); for all frequencies (s.f)
            %[Y11, Y12, Y21, Y22] = toYparams(freq); for a specified frequency
            if obj.numPorts ~= 2
				warning('Only for 2 port measurements');
            end
            switch nargin
				case 2
					fS11 = obj.S11; fS12 = obj.S12; fS21 = obj.S21; fS22 = obj.S22;
				case 3
					S = obj.toMat(freq);
					fS11 = S(1, 1); fS12 = S(1, 2); fS21 = S(2, 1); fS22 = S(2, 2);
            end
            delS = (1 + fS11).*(1 + fS22) - fS12.*fS21;
            Y11 = ((1 - fS11).*(1 + fS22) + fS12.*fS21)./delS./obj.Z0;
            Y12 = -2.*fS12./delS./obj.Z0;
            Y21 = -2.*fS21./delS./obj.Z0;
            Y22 = ((1 + fS11).*(1 - fS22) + fS12.*fS21)./delS./obj.Z0;
        end

		function writeCSV(obj, filename)
            %Write a CSV file in the format
            %f, S11 - for 1 port
            %f, S11, S21, S12, S22 - for 2 port
            %f, S11, S12, S13, S21, S22, S23, S31, S32, S33 - for 3 port
            %And so on
			switch nargin
				case 1
					[~, csvname, ~] = fileparts(obj.filename);
					csvwrite(strcat(csvname, '.csv'), obj.data);
				case 2
					csvwrite(filename, obj.data);
			end
        end
        
        %Work in progress
%         function writeSNP(obj, filenameOut)
%             fileExt = ['.s' num2str(obj.numPorts), 'p'];
%             
% 			dlmwrite(filenameOut, obj.data, ' ');
%         end
        
        function cascade(obj, N)
            %Calculates S parameters when cascaded N times. Only valid for
            %2 port unit cells.
            %This function will overwrite the original S-parameters, so it
            %may be advisable to copy the original unit cell first using
            %copyobj()
                [Au, Bu, Cu, Du] = obj.toABCD;
                At = Au; Bt = Bu; Ct = Cu; Dt = Du;
                for i = 1:N-1
                    Attemp = At; Bttemp = Bt; Cttemp = Ct; Dttemp = Dt;
                    At = Attemp.*Au + Bttemp.*Cu;
                    Bt = Attemp.*Bu + Bttemp.*Du;
                    Ct = Cttemp.*Au + Dttemp.*Cu;
                    Dt = Cttemp.*Bu + Dttemp.*Du;
                end
                obj.S11 = (At + Bt/obj.Z0 - Ct.*obj.Z0 - Dt)./(At + Bt/obj.Z0 + Ct.*obj.Z0 + Dt);
                obj.S12 = (2.*(At.*Dt - Bt.*Ct))./(At + Bt/obj.Z0 + Ct.*obj.Z0 + Dt);
                obj.S21 = 2./(At + Bt/obj.Z0 + Ct.*obj.Z0 + Dt);
                obj.S22 = (-At + Bt./obj.Z0 - Ct.*obj.Z0 + Dt)./(At + Bt/obj.Z0 + Ct.*obj.Z0 + Dt);
        end
        
        function plotDispersion(obj, flip)
            %Plots the dispersion diagram when s is a unit cell.
            %s.plotDispersion();
           Beta_p = real(acosd((1 - obj.S11.*obj.S22 + obj.S12.*obj.S21)./(2.*obj.S21)));
           %Beta_p = acos((1 - S(1,1)*S(2,2) + S(1,2)*S(2,1))/(2*S(2,1)))
           switch nargin
               case 1
                    plot(Beta_p, obj.f/obj.freqScale{1}, 'k');
                    ylabel(['f (' obj.freqScale{2} ')'])
                    xlabel('\beta p (deg)')
               case 2
                   if ~strcmpi(flip, 'flip')
                       plot(Beta_p, obj.f/obj.freqScale{1}, 'k');
                       ylabel(['f (' obj.freqScale{2} ')'])
                       xlabel('\beta p (deg)')
                   else
                       plot(obj.f/obj.freqScale{1}, Beta_p, 'k');
                       xlabel(['f (' obj.freqScale{2} ')'])
                       ylabel('\beta p (deg)')
                   end
           end

           grid on
        end %
        

    end

    methods(Static)
        
        function sDUT = deembed(err1, sALL, err2)
            %Deembed S-parameters under the following conditions
            %Measurement: sALL = P1 <-- err1 <--> DUT <--> err2 --> P2
            %sDUT = SPARAMS.deembed(err1, sALL, err2);
            %Will return a new sparameter object that represents the DUT
            %CURRENTLY UNDER TESTINNG
            if err1.f ~= err2.f | err1.f ~= sALL.f
                error('frequencies unequal');
            end
            len = length(sALL.f);
            
            Tt11 = zeros(1, len); Tt12 = zeros(1, len); Tt21 = zeros(1, len); Tt22 = zeros(1, len); 
           for i = 1:len
               f = sALL.f(i);
               [e1T11, e1T12, e1T21, e1T22] = err1.toTParams(f);
               [dT11, dT12, dT21, dT22] = sALL.toTParams(f);
               [e2T11, e2T12, e2T21, e2T22] = err2.toTParams(f);
               tempT = [e1T11 e1T12;e1T21 e1T22]\[dT11 dT12;dT21 dT22]*inv([e2T11 e2T12;e2T21 e2T22]);
               Tt11(i) = tempT(1, 1);
               Tt12(i) = tempT(1, 2);
               Tt21(i) = tempT(2, 1);
               Tt22(i) = tempT(2, 2);               
           end
           
           warning('off', 'all')
           sDUT = SPARAMS();
           sDUT.setNumPorts(2);
           warning('on', 'all')
           sDUT.S11 = Tt12./Tt22;
           sDUT.S12 = (Tt11.*Tt22 - Tt12.*Tt21)./Tt22;
           sDUT.S21 = 1./Tt22;
           sDUT.S22 = -Tt21./Tt22;
           sDUT.f = sALL.f;
           sDUT.setZ0(sALL.Z0);

        end
        
        function plotAnydB(ob1, splt1, ob2, splt2)
            %Plot any S-parameters from any 2 objects in dB
            %SPARAMS.plotAnydB(obj1, {'S11', 'S21'}, obj2, {'S11', 'S21'});
            if class(splt1) == 'char' | class(splt2) == 'char'
                error('Input name of parameters in cell format')
            end

            figure();
            hold on;
            splt1 = cell(splt1); splt2 = cell(splt2);
            for i = 1:length(splt1)
                stringToPlot = ['ob1.' splt1{i}];
                stoPlot = eval(stringToPlot);
                plot(ob1.f, 20*log10(abs(stoPlot)));
            end
            for i = 1:length(splt2)
                stringToPlot = ['ob2.' splt2{i}];
                stoPlot = eval(stringToPlot);
                plot(ob2.f, 20*log10(abs(stoPlot)));
            end
            xlabel('f (Hz)');
            ylabel('S-parameters (dB)');
        end

		function Z = s2z(S, Z0)
            %General conversion from S to Z
			fS11 = S(1, 1); fS12 = S(1, 2); fS21 = S(2, 1); fS22 = S(2, 2);
			Z11 = Z0*((1 + fS11).*(1 - fS22) + fS12.*fS21)./((1 - fS11).*(1 - fS22) - fS12.*fS21);
			Z12 = Z0*(2.*fS12)./((1 - fS11).*(1 - fS22) - fS12.*fS21);
			Z21 = Z0*(2*fS21)./((1 - fS11).*(1 - fS22) - fS12.*fS21);
			Z22 = Z0*((1 - fS11).*(1 + fS22) + fS12.*fS21)./((1 - fS11).*(1 - fS22) - fS12.*fS21);
			Z = [Z11 Z12;Z21 Z22];
		end

		function S = z2s(Z, Z0)
            %General conversion from Z to S
			Z11 = Z(1, 1); Z12 = Z(1, 2); Z21 = Z(2, 1); Z22 = Z(2, 2);
			dZ = (Z11 + Z0).*(Z22 + Z0) - Z12.*Z21;
			S11 = ((Z11 - Z0).*(Z22 + Z0) - Z12.*Z21)/dZ;
			S12 = 2*Z12.*Z0./dZ;
			S21 = 2*Z21.*Z0./dZ;
			S22 = ((Z11 + Z0).*(Z22 - Z0) - Z12.*Z21)/dZ;s
			S = [S11 S12;S21 S22];
		end
    end
end

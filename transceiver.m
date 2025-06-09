classdef transceiver < handle
    properties
        b
        ss
        sc
        sp
        t
        fc
        Fs
        T
        levels
        map
        modulation_type
        num_bits
    end

    methods
        function data = load_data_from_file(~, file_path)
            % INPUT : path to data file.
            % the file must have binary values seperated by spaces
            % OUTPUT : binary data as an array

            file = fopen(file_path, 'r');   
            if file ==-1   % file does not exist
                disp('Invalid file path')
            end
            data = fscanf(file, '%d');  % read data
            data = data';               % transform column vector to row vector
            fclose(file);
        end
        
        %% MODULATION SCHEMES
        
        % INPUT : data and no. of bits per symbol
        % OUTPUT: encoded data array

        
        % 1. PAM modulation
        function [symbols, symbols_imaginary] = PAM_modulation(self, data, b)
           M = 2^b;  % number of levels
           self.levels = -(M-1): 2: M-1;  
           self.map = 0:1:M-1;  % 0 -> -(M-1), 1-> -(M-2), ...., M-1 -> (M-1)

           p = 1; k = 0;
           symbols = [];
           for i= 1:length(data)
                k = k + p*data(i);
                p = p*2;
                % after every b bits, assigns symbol for the bits
                if mod(i, b)==0 | i ==length(data) 
                    symbols = [symbols, k];
                    k = 0; p = 1;
                end
           end
           
           symbols = self.levels(1 + symbols);
           symbols_imaginary = zeros(1, length(symbols));  % imaginary part is always zero
        end
        
        % 2. QAM modulation
        function [symbols_real, symbols_imaginary] = QAM_modulation(self, data, b)
            M = 2^(b/2);  % Number of levels per dimension (sqrt of constellation size)
            self.levels = -(M-1):2:(M-1);  % PAM levels, e.g., [-3 -1 1 3] for M=4
            symbols_real = [];    % first half bits are mapped to the real part
            symbols_imaginary = [];  % second half bits are mapped to the imaginary part
            
            k = 0; p = 1;
            for i = 1:length(data)
                k = k + p * data(i);  % Accumulate bits into integer
                p = p * 2;
        
                if mod(i, b) == 0 || i == length(data)
                    a = mod(k, M);                 % I-component index
                    bval = floor(k / M);           % Q-component index
                    symbols_real = [symbols_real, self.levels(a + 1)];
                    symbols_imaginary = [symbols_imaginary, self.levels(bval + 1)];
                    k = 0; p = 1;  % Reset for next symbol
                end
            end
        end
        
        % 3. QPSK modulation
        function [symbols_real, symbols_imaginary] = QPSK_modulation(self, data, b)
            % QPSK = 2 bit QAM
            [symbols_real, symbols_imaginary] =  self.QAM_modulation( data, 2); 
        end
        
        %% Decoders
        % Function : given an modulation value, gives the corresponding
        % binary value
        % INPUT : real and imaginary parts, no. bits used for encoding
        % OUTPUT: integer value of the binary sequence

        % 1. QAM decoder
        function z = QAM_decoder(self, x, y, b)
            M = 2^(b/2);
            levels = -(M-1):2:(M-1);
            [~, idx_x] = min(abs(x - levels)); 
            [~, idx_y] = min(abs(y - levels));
            z = (idx_y-1)*M + idx_x-1;
        end
        % 2. PAM decoder
        function z = PAM_decoder(self, x, y, b)
            M = 2^(b);
            levels = -(M-1):2:(M-1);
            [~, idx_x] = min(abs(x - levels));
            z = idx_x-1;
        end
        % 3. QPSK decoder
        function z = QPSK_decoder(self, x, y, b)
            z = self.QAM_decoder(x, y, 2);
        end
        
        function display_modulation(self, modulation_type, b)
            % FUNCTION : displays modulation constellation
            % INPUT : modulation type, number of bits
            % OUTPUT : plot
            if strcmp(modulation_type, "PAM")
                encoder = @(x, y) self.PAM_modulation(x, y);
                decoder = @(x, y, b) self.PAM_decoder(x, y, b);
            elseif strcmp(modulation_type, "QAM")
                encoder = @(x, y) self.QAM_modulation(x, y);
                decoder = @(x, y, b) self.QAM_decoder(x, y, b);
            else
                encoder = @(x, y) self.QPSK_modulation(x, y);
                decoder = @(x, y, b) self.QPSK_decoder(x, y, b);
                b = 2;
            end
            N = 2^b;  % Total number of combinations
            bit_array = zeros(1, N * b);  % Preallocate 1D array
            
            % generates array all values that can be represented by b bits
            for i = 0:N-1
                bits = bitget(i, 1:b);              % LSB to MSB
                bit_array(i*b + 1 : (i+1)*b) = bits;
            end
            
            [bx, by] = encoder(bit_array, b);  % apply encoding
            
            % plot properties
            scatter(bx, by, 100,"rx");  %scatter plot
            xlabel('In-phase');  ylabel('Quadrature');
            yticks(min(by)-2:2:max(by)+2); xticks(min(bx)-2:2: max(bx)+2)
            grid on; axis equal;
            xline(0, '-k'); yline(0, '-k'); xlim([min(bx)-2 max(bx)+2]); ylim([min(by)-2 max(by)+2])
            for i = 1:length(bx)
                label = decoder(bx(i), by(i), b);  % binary value corresponding to each point
                text(bx(i)-0.2, by(i) + 0.4, dec2bin(label,b), 'FontSize', 10);
            end
            xlabel("Real axis","interpreter","latex"); ylabel("Complex axis","interpreter","latex"); 
            title(modulation_type + " constellation");
        end
        
        %% pulse shaping filter
        function p = p_rect_NRZ(~, t, T)
            % rectangle function that stays on for time 0 to T
            p = double(t<T & t>0);
        end
        
        %% PULSE SHAPING AND BANDPASS CONVERSION
        function sp = modulation(self, data, modulation_type, num_bits, Fs, T, fc)
            % INPUT : binary data, modulation type, 
            %         no. of bits per symbol, band pass frequency, 
            %         samp. frequency, Bit rate, bandpass frequency
            % OUTPUT : Transmitted signal
            self.modulation_type = modulation_type;
            self.num_bits = num_bits;
            

            if strcmp(modulation_type, "PAM")   % chose modulation function
                encoder = @(x, y) self.PAM_modulation(x, y);
            elseif strcmp(modulation_type, "QAM")
                encoder = @(x, y) self.QAM_modulation(x, y);
            else
                encoder = @(x, y) self.QPSK_modulation(x, y);
                self.num_bits = 2;
            end
            
            % find the symbol values (real and imaginary parts)
            [bx, by] = encoder(data, num_bits);
            self.b = [bx; by];
                
            Tmin = -5*T;
            Tmax = length(self.b)*T + 5*T;
            t = Tmin : (1/Fs) : Tmax;  %  time axis
            
            self.Fs = Fs; self.t = t; self.T = T;
            self.ss = zeros(1, length(t));
            self.sc = zeros(1, length(t));

            i = 1; t_cur = 0;
            
            % sc: in phase component,  ss: quadrature component
            while i<=length(bx)
                % The amplitudes are multiplied with the shaping function
                % over time intervals of size T
                self.sc = self.sc + bx(i)* self.p_rect_NRZ(t-(i-1)*T, T);
                self.ss = self.ss + by(i)* self.p_rect_NRZ(t-(i-1)*T, T);
                t_cur = t_cur+T;
                i = i+1;
            end
            
            % Baseband to Passband conversion

            self.fc = fc;
            I = self.sc.*cos(2*pi*fc*t);  % in phase component 
            Q = self.ss.*sin(2*pi*fc*t);  % Quadrature component

            sp = sqrt(2)*(I-Q);  % passband channel waveform
            self.sp = sp;
            
        end
        
        %% DEMODULATION FILTERS
        % INPUT: received signal
        % OUTPUT: in-phase and quadrature symbol amplitudes

        % 1. Lowpass filter
        function r = lowpass_filter(self, rp)
            t = self.t; fc = self.fc; % synchronizing with transmitter's parameters
            Fs = self.Fs; T = self.T;
            
            % multiply with cosine and sine function of pass band frequency
            % and apply low pass filter to obtain the baseband signals
            rI = (sqrt(2))*rp.*cos(2*pi*fc*t); 
            rQ = -(2^0.5)*rp.*sin(2*pi*fc*t);
            rc = lowpass(rI, fc/10, Fs);
            rs = lowpass(rQ, fc/10, Fs);
            [f, Rp] = self.calc_fft(rp, Fs);
            [f, Rc] = self.calc_fft(rc, Fs);
            subplot(2, 1, 1)
            plot(t, rc); xlabel("time","interpreter","latex"); ylabel("rc","interpreter","latex"); title("Filtered in-phase component"); 
            subplot(2, 1, 2)
            plot(f, abs(Rc)); xlabel("frequency","interpreter","latex"); ylabel("Rc","interpreter","latex"); title("FFT of filtered in-phase component");
            
            
            % sample the waveform at time interval T to obtain the amplitudes
            ind = round((T/5 - t(1))*Fs);
            t_cur = T/5;
            N = round(t(end)/T);
            r = zeros(2,N);
            i = 1;
            while t_cur< t(end)
                r(1, i) = rc(ind);
                r(2, i) = rs(ind);
                t_cur = T + t_cur;
                i = i+1;
                ind = round(ind + T*Fs);
            end
        end
        
        % 2. Correlator
        function r = correlator(self, rp)
            fc = self.fc; Fs = self.Fs;
            T = self.T; t = self.t;

            N = round(t(end)/T);
            r = zeros(2,N);
            t_cur = 0; i = 1;
            
            % calculate inner product with the orthogonal components -
            % p(t)sin(wt) and p(t)cos(wt) to obtain the coefficients
            while t_cur < t(end)
                p = self.p_rect_NRZ(t-(i-1)*T, T);
                x = sum(rp.*p.*cos(2*pi*t*fc))*(sqrt(2)/Fs); % in-phase
                y = sum(-rp.*p.*sin(2*pi*t*fc))*(sqrt(2)/Fs); % quadrature
                r(1,i) = x;
                r(2,i) = y;
                i=i+1;
                t_cur = t_cur+T; % shift time by T
            end
        end
        
        % 3. Matched Filter
        function r = matched_filter(self, rp)
            fc = self.fc; Fs = self.Fs;
            T = self.T; t = self.t;
            
            N = round(t(end)/T); % number of symbols
            r = zeros(2, N);
            
            % Define matched filter (time-reversed pulse shape)
            p = self.p_rect_NRZ(t, T);   % assuming t is symmetric around 0
            mf = fliplr(p);              % matched filter = time-reversed pulse
        
            % I and Q baseband signals
            rI = sqrt(2)*rp .* cos(2*pi*fc*t);
            rQ = -sqrt(2)*rp .* sin(2*pi*fc*t);
        
            % Matched filtering
            x = conv(rI, mf, 'full') / Fs;
            y = conv(rQ, mf, 'full') / Fs;
            plot(x)
            xlabel("time","interpreter","latex"); ylabel("signal","interpreter","latex"); title("Matched filter output signal");
            self.fix_plot_dim(250);

            % Sample at symbol centers
            sample_indices = round( (0:N)*T*Fs + (t(end)-t(1))*Fs + 1 ); % assuming t(1) ~ 0
            for i = 1:N
                if sample_indices(i) <= length(x)
                    r(1, i) = x(sample_indices(i));
                    r(2, i) = y(sample_indices(i));
                else
                    r(:, i) = 0;
                end
            end
        end
        
        % adds white guassian noise to the channel at given SNR (db)
        function sp = add_noise(self, sp, SNR)
           sp = awgn(sp, SNR, 'measured');
        end
        
        %% Demodulation
        function r = demodulation(self, rp, type)
            % INPUT : received signal, demodulation type
            % output : symbol amplitude array
            
            if strcmp(type, 'lowpass')
                demod_func = @(x) self.lowpass_filter(x);
            elseif strcmp(type, 'correlator')
                demod_func = @(x) self.correlator(x);
            else 
                demod_func = @(x) self.matched_filter(x);
            end

            r = demod_func(rp);
        end
        
        
        function d = get_binary_sequence(self, r)
            % INPUT - symbol amplitudes
            % OUTPUT - decoded binary array

            % PAM
            if strcmp(self.modulation_type, "PAM")
                 N = size(r, 2);
                 c = zeros(1, N);
                 for i = 1:N;
                    [~, idx] = min(abs(r(1, i) - self.levels));
                    c(i) = self.map(idx);
                 end
                 d = [];
                 for i = 1:N
                     bin_str = dec2bin(c(i), self.num_bits);
                     binary_array = bin_str - '0';
                     d = [d, flip(binary_array)];
                 end

            % QAM
            elseif strcmp(self.modulation_type, "QAM")
                d = [];
                N = size(r, 2);
                
                for i = 1:N;
                    z = QAM_decoder(self, r(1,i), r(2,i), self.num_bits);
                    bin_str = dec2bin(z, self.num_bits);       
                    bit_array = bin_str - '0';
                    d = [d, flip(bit_array)];
                end
            
            % QPSK
            else
                d = [];
                N = size(r, 2);
                
                for i = 1:N;
                    z = QPSK_decoder(self, r(1,i), r(2,i), 2);
                    bin_str = dec2bin(z, self.num_bits);       
                    bit_array = bin_str - '0';
                    d = [d, flip(bit_array)];
                end
            end
        end
        
        function demodulation_constellation(self, r)
            scatter(r(1,:), r(2,:),200,'b.');
            hold on;
            self.display_modulation(self.modulation_type, self.num_bits);
            % hold on;
            title(self.modulation_type + " demodulation constellation")
            N = size(r, 2);
        end
        function fix_plot_dim(~, w)
            fig = gcf;
            current_position = fig.Position;
            fig.Position = [current_position(1), current_position(2), current_position(3), w];
        end

        function plot_data(~,t, x,XLABEL,YLABEL, TITLE) 
            plot(t,x)           % plot data with legend
            xlabel(XLABEL,"interpreter","latex")      % label with latex font
            ylabel(YLABEL,"interpreter","latex")                      % update legend
            title(TITLE)                              % Set title
            grid on;                                   % draw grid
        
            %fix height of plot for comfortable viewing
            fig = gcf;
            current_position = fig.Position;
            fig.Position = [current_position(1), current_position(2), current_position(3), 250];
        end

        function [f, X] = calc_fft(~,x, Fs)
            n = length(x);
            X = fft(x);
            f = (0:n-1)*(Fs/n);
        end
    end
end
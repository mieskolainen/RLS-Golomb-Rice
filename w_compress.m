% RLS-Rice/Golomb lossless encoder (fixed to 16 bit integers)
%
% w_compress(input_wave, output_binary, M, b)
%
% INPUT:     input_wave  :  Wave file path, e.g. '/home/test.wav'
%         output_binary  :  Output binary path, e.g. '/home/codec.bin' 
%                     M  :  Prediction order for RLS, (1...32)
%                     b  :  Block size, e.g. 1024, in powers of two            
%
%
% Mikael Mieskolainen, 2011

function w_compress(input_wave, output_binary, M, b)

% READ THE WAVE FILE

[y_orig, Fs] = audioread(input_wave);

% Scale the signal to 16-bit integer range.
y = y_orig * 2^15;


% RLS (RECURSIVE LEAST SQUARES)

% Initializations of RLS-algorithm.
L = length(y);
e = zeros(1,L);
P = eye(M);
w = zeros(M,1);
lambda = 0.999;

% RLS decorrelation
for n = M+1:L
    
    if (mod(n,10000) == 0)
        fprintf('RLS encoding on sample: %d / %d \n', n, L);
    end
    
    % Select a short vector of signal y, denote this as y_v
    y_v = y(n-1:-1:n-M);
    [P, w, e(n)] = RLSo(y_v, P, w, y(n), lambda);
end


% RICE-GOLOMB BINARY ENCODING

% Open the file
file = fopen(output_binary, 'w');

% Remove the first M errors from e (because they are zero)
e = e(M+1:end);

% Make sure we have a signal length multiple of b (if necessary)
if (mod(length(e), b) ~= 0)
    e = [e zeros(1, b - mod(length(e),b))];
end

% Write the header information [e_length, y_length, M, b, y(M:-1:1), Fs]
fwrite(file, length(e), 'uint32');
fwrite(file, length(y), 'uint32');
fwrite(file, M, 'uint8');
fwrite(file, b, 'uint16');
fwrite(file, y(M:-1:1), 'int16');
fwrite(file, Fs, 'uint32');


% ENCODING LOOP
for i = 1:b:length(e)
    
    if (mod(i,10) || mod(i,21))
        fprintf('Golomb-Rice encoding on sample: %d / %d \n', i, length(e));
    end
    
    % This data block
    data = e(i:i+b-1);
    
    % Check all the different p's (0 ... 15) in a brute force style
    % p means, how many significant bits
    wl = zeros(16,1);
    for p = 0:15
        
        % Go through all the values x in data
        wls = 0; % Word lengths
        for j = 1:length(data)
            wls = wls + (1 + p + (floor(abs(data(j))/(2^p)) + 1));
        end
        wl(p+1) = wls;
    end
    
    % Find the best p
    [~, p] = min(wl);
    p = p - 1;                % Fix the Matlab indexing
    
    % Code the p to file
    fwrite(file, p, 'ubit4'); % Write the parameter p value for this block
    
    % Code the values from this block with the best p
    for j = 1:length(data)
        colombenc(file, data(j), p);
    end
end

fclose(file);

end


% Recursive Least Squares function for encoding
% ------------------------------------------------------------------------
%
% Input:      y       =  Input vector [y(n-1), y(n-2), ..., y(n-M)], (Mx1)
%             P       =  Matrix P_{n-1},         (MxM)
%             w       =  Weight vector w_{n-1},  (Mx1)
%             y_n     =  Sample y(n),            (1x1)
%             lambda  =  Forgetting factor,      (1x1)
%
% Output:     P       =  Matrix P_{n},           (MxM)
%             w       =  Weight vector w_{n},    (Mx1)
%             e       =  Prediction error,       (1x1)
%
% Mikael Mieskolainen, 2011

function [P, w, e] = RLSo(y, P, w, y_n, lambda)

% Vector k
k = (lambda^(-1) * P * y) / (1 + lambda^(-1) * y' * P * y);

% Estimate signal value (scalar)
y_hat = round(w'*y);

% Error value (scalar)
e = y_n - y_hat;

% New weight vector
w = w + e*k;

% New P
P = lambda^(-1) * P - lambda^(-1) * k * y' * P;

end


% Golomb-Rice code encoder
% ------------------------------------------------------------------------
% 
% Input:    file  =  Filehandle
%              S  =  Value to be coded
%              p  =  Parameter p
%
% 
% Any symbol S can be represented as a quotient (Q) and remainder (R),
% where S = Q x M + R
%
% and M = 2^p  <-> log2 (M) = p
%
% Mikael Mieskolainen, 2011

function colombenc(file, S, p)

sign_bit = heavi(S);        % Take the sign bit (0,1)
S = abs(S);                 % Take the absolute value

R = bitand(S, 2^p - 1);     % The least significant bits (LSB)
Q = floor(S/(2^p));         % The most significant bits (MSB)

% Now write all the bits
fwrite(file, sign_bit, 'ubit1');  % Write sign bit (0,1)

if (p ~= 0)
    fwrite(file, R, ['ubit' int2str(p)]); % LSB bits (remainder)
end

% Write MSB bits (quotient) using unary coding
n = Q;
while (n > 0)
    fwrite(file, 1, 'ubit1');     % Write bit 1
    n = n - 1;
end
fwrite(file, 0, 'ubit1');         % Write bit 0

end

% Heaviside function
function out = heavi(in)

if (in >= 0)
    out = 1;
else
    out = 0;
end

end
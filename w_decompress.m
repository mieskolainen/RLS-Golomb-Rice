% RLS-Rice/Golomb lossless decoder for 16 bit integers
%
% w_decompress(input_binary, output_wave) 
%
% INPUT:    input_binary  :  Input binary file, e.g '/home/test.bin'
%            output_wave  :  Output wave file, e.g '/home/output.wav'
%
%
% Mikael Mieskolainen, 2011

function w_decompress(input_binary, output_wave)


% GOLOMB-RICE DECODING

% Open the file
file = fopen(input_binary, 'r');

% Read the header information: [e_length, y_length, M, b, y(M:-1:1), Fs]
e_length = fread(file, 1, 'uint32');
y_length = fread(file, 1, 'uint32');
M = fread(file, 1, 'uint8');
b = fread(file, 1, 'uint16');
y_v = fread(file, M, 'int16');
Fs = fread(file, 1, 'uint32');

% We have M zeros at the beginning of e, so add M
e_length = e_length + M;


% Create the error signal vector
e = zeros(1, e_length);

% Decode Golomb-Rice, start after M zeros of the beginning of e
k = 1;
for i = M+1:b:e_length
    
    if (mod(k,100) == 0)
        fprintf('Colomb-Rice decoding on sample %d / %d \n', i, e_length);
    end
    e(i:i+b-1) = golombdec(file, b)';
    k = k + 1;
end

% Close the file
fclose(file);


% RLS DECODER

% Decoded final signal vector
yd = zeros(e_length, 1);
yd(1:M) = y_v(M:-1:1);

P = eye(M);
w = zeros(M,1);
lambda = 0.999;

for n = M+1:length(yd)
    
    if (mod(n,10000) == 0)
        fprintf('RLS decoding on sample: %d / %d \n', n, length(yd));
    end
    
    if (n ~= M+1)
       y_v(2:end) = y_v(1:end-1);  % Shift the values for the new one
       y_v(1) = yd(n-1);           % Select the previous decoded y(n)
    end
    
    % Run the RLS
    [P, w, y_hat] = RLSd(y_v, P, w, e(n), lambda);
    
    % Scalar value y(n) = e(n) + y_hat(n)
    yd(n) = e(n) + y_hat;
end

% Final original y
y = yd(1:y_length);

% WRITE THE WAVE FILE
audiowrite(output_wave, int16(y), Fs, 'BitsPerSample', 16);

end


% Golomb-Rice code decoder
% ------------------------------------------------------------------------
% 
% Input:    file  =  Filehandle
%              b  =  Block size (length)
%
%
% Any symbol S can be represented as a quotient (Q) and remainder (R),
% where S = Q x M + R
%
% and M = 2^p  <-> log2 (M) = p
%
% Output:   output = decoded block
%
% Mikael Mieskolainen, 2011

function output = golombdec(file, b)

p = fread(file, 1, 'ubit4');            % Read the parameter p for this block
output = zeros(b,1);                    % Output vector (block values)

for i = 1:b
    
    sign_bit = fread(file, 1, 'ubit1');         % Read sign bit (0,1)
    if (sign_bit == 0), sign_value = -1; else sign_value = 1; end
    
    if (p ~= 0)
       R = fread(file, 1, ['ubit' int2str(p)]); % Read the LSB bits
    else
       R = 0;
    end
    
    % Unary decode the quotient Q from file
    Q = 0;
    while (fread(file,1,'ubit1') == 1)          % Read while we get ones
        Q = Q + 1;
    end
    
    % Write the final value
    output(i) = sign_value * (Q * (2^p) + R);
end

end


% Recursive Least Squares function for decoding
% ------------------------------------------------------------------------
%
% Input:      y       =  Input vector [y(n-1), y(n-2), ..., y(n-M)], (Mx1)
%             P       =  Matrix P_{n-1},         (MxM)
%             w       =  Weight vector w_{n-1},  (Mx1)
%             e       =  Error value e(n)        (1x1)
%             lambda  =  Forgetting factor,      (1x1)
%
% Output:     P       =  Matrix P_{n},           (MxM)
%             w       =  Weight vector w_{n},    (Mx1)
%             y_hat   =  Predicted value,        (1x1)
%
% Mikael Mieskolainen, 2011

function [P, w, y_hat] = RLSd(y, P, w, e, lambda)

% Vector k
k = (lambda^(-1) * P * y) / (1 + lambda^(-1) * y' * P * y);

% Estimate signal value (scalar)
y_hat = round(w'*y);

% New weight vector
w = w + e*k;

% New P
P = lambda^(-1) * P - lambda^(-1) * k * y' * P;

end
% Lossless compression based on Recursive Least Squares and Colomb-Rice
% coding
%
% Mikael Mieskolainen, 2011
clear; close all;


%% Part I

[y_orig, Fs] = audioread('mike.wav');
fprintf('\nTask 1:\n');

% Scale the signal to 16-bit signed integer range.
y = y_orig*2^15;

% Normalization of histogram samples. Sample space is 2^16, because that's
% the cardinality of the discrete levels of signal y.
[p, xout] = hist(y, -2^15:2^15-1);
p = p./sum(p);

% Plot of the normalized histogram.
bar(xout, p); set(gca,'yscale','log'); axis([-2^14 2^13 0 inf]);
xlabel('Symbol (Amplitude)'); ylabel('Probability p(x)');

% Calculate the entropy.
p = p(p ~= 0); % Take only non-zero bins
H = -sum(p .* log2(p));
fprintf('Empirical entropy for signal y is %0.2f.\n', H)
%}


%% Part II

fprintf('\nPart 2:\n');

% Use 7zip for compression. The command used for compression:
% 7z a -mx9 mike.7z mike.wav
d1 = dir('mike.7z');
d2 = dir('mike.wav');
compression_ratio = d1.bytes / d2.bytes;
fprintf('Compression ratio with 7zip is %0.3f.\n', compression_ratio)


%% Part III

M = 16; %[4 8 16 32 64]; % Prediction order
b = 32; %[32 64 128 256 512 1024]; % Block size

for i = 1:length(M)
    for j = 1:length(b)
        % Compress file
        w_compress('mike.wav', sprintf('compressed_M%d_b%d.bin', M(i), b(j)), M(i), b(j));
    end
end


%% Part IV

% Decompress file and save as wav
w_decompress('compressed_M16_b32.bin', 'decompressed.wav');

% Read out the original and compressed/decompressed version
[orig, Fs] = audioread('mike.wav');
[comp, ~]  = audioread('decompressed.wav');

% Check that they are exactly same => Lossless
diff = norm(orig - comp)

d1 = dir('testfile_M16_b32.bin');
d2 = dir('mike.wav');

comp_ratio = d1.bytes / d2.bytes;
fprintf('Compression ratio with the RLS codec is %0.3f \n', comp_ratio);
%}


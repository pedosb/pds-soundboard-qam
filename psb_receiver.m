function [err, snr, preErr, wOffset, pOffset, rxSignalBlock] = psb_receiver(...
    folder, trStart, ignoreN, pBlocks, txGain, debug)
%% psb_receiver Decode and evaluate a PAM transmission in blocks
% Input variables:
%   folder:  Folder with data.mat (decoding variables), data.bin
%            (transmitted signal) and rx.bin (received signal)
%   trStart: Where to search for the first block
%   pBlocks: Process pBlocks blocks, default is all blocks
%   ignoreN: Does not consider the last ignoreN symbols from payload
% Output variables:
%   err:     Number of different symbols on payload
%   snr:     SNR for payload
%   preErr:  Number of different symbols on preamble
%   wOffset: Estimated frequency offset
%   pOffset: Estimated phase offset

if nargin < 6
    debug = 0;
end
if nargin < 5
    txGainStr = '';
else
    txGainStr = sprintf('_%.2d', txGain);
end
if nargin < 4
    pBlocks = 0;
end

rxFN = [folder 'rx' txGainStr '.bin'];
txIndicesFn = [folder 'tx_indices.bin'];
txSignalFn = [folder 'data.bin'];

% Load decoding variables and rename some of them
load([folder 'data.mat']);
L = oversamplingFactor;
M = modulationOrder;
gDelay = 50;
preambleLength = length(txSignalPreamble) / L;
if pBlocks ~= 0
    nBlocks = pBlocks;
end

nTrG = 1;
% Number of symbol error on each block
err = zeros(nBlocks, payloadLenght - ignoreN, nTrG);
snr = zeros(nBlocks, 1, nTrG);     % SNR of each block
preErr = zeros(nBlocks, 1, nTrG);  % Preamble error of each block
wOffset = zeros(nBlocks, 1, nTrG); % Frequency offset of each block
pOffset = zeros(nBlocks, 1, nTrG); % Phase offset of each block
% Received payload signal without frequency and phase offset correction
rxSignalBlock = zeros(nBlocks, payloadLenght, nTrG);
for trG=1:nTrG
    if trG ~= 1
        trStart = startOfFile + preStart +... % Use last prediction
            guardUsrpSamplesLength + 1 * ((blockLenght * L) / 2); % Go to half a block
    end
for bn=1:nBlocks;
    % Skip some samples (analog USRP stabilization)
    startOfFile = trStart + (bn - 1) * blockLenght * L;

    %% Load TX data
    % Load transmitted block considering the tx filter group delay and the
    % guard samples
    txSignal = read_complex(txSignalFn, blockLenght*L,...
        guardUsrpSamplesLength+gDelay+(bn-1)*blockLenght*L, 'float32');
    % Sample symbols from transmitted signal
    txSymbols = txSignal(preambleLength*L+1:L:end);
    % Load transmitted indices for the current block
    txIndices = read_array(txIndicesFn, payloadLenght,...
        (bn-1)*payloadLenght, 'uint8');

    %% Find the start of a block
    rxSignal = read_complex(rxFN, blockLenght*L, startOfFile);
    
    rxSignal = correct_offset(rxSignal, 1:length(rxSignal));
    
    [preCorr, preCorrLags] = xcorr(rxSignal, txSignalPreamble);
    % TODO: Warn if preCorrMax is too small
    % Find max at cross correlation
    [preCorrMax, preCorrMaxIdx] = max(preCorr);
    % Find index at rxSignal where the preamble starts
    preStart = preCorrLags(preCorrMaxIdx);

    %% Process a block
    % Read the block
    rxBlock = read_complex(rxFN, blockLenght*L, startOfFile+preStart);
    rxSignalBlock(bn, :, trG) = rxBlock(preambleLength*L+1:L:end);
    [rxBlock, wOffset(bn, 1, trG), pOffset(bn, 1, trG)] =...
        correct_offset(rxBlock, 1:preambleLength*L);
    
    % Downsample block to symbol rate
    rxBlockSymbols = rxBlock(1:L:end);
    % Extract preamble
    rxPreambleSymbols = rxBlockSymbols(1:preambleLength);
    % Extract payload
    rxSymbols = rxBlockSymbols(preambleLength+1:end);

    %% Correct pi ambiguity
    % The phase offset estimated by FFT has an ambiguity in pi
    rot = 0;
    if real(sum(txPreambleSymbols .* rxPreambleSymbols)) < 0
        rot = pi;
    end

    rxPreambleSymbols = rxPreambleSymbols .* exp(-1j*rot);
    rxSymbols = rxSymbols .* exp(-1j*rot);
    pOffset(bn, 1, trG) = pOffset(bn, 1, trG) + rot;

    %% Normalize symbols power
    % Power is 1 for a 2-PAM preamble
    rxPreambleSymbols =...
        (sqrt(1) / rms(rxPreambleSymbols)) * rxPreambleSymbols;
    % Normalize payload energy in nNormBlock symbols blocks
    constRms = rms(qammod(0:M-1, M));
    txSymbols = (constRms / rms(txSymbols)) * txSymbols;
    % Normalize block to the theoretical PAM power
    rxGain = (constRms / rms(rxSymbols));
    rxSymbols = rxGain * rxSymbols;
    rxSignalBlock(bn, :, trG) = rxGain * rxSignalBlock(bn, :, trG);

    %% Performance evaluation
    disp(['Block ' num2str(bn)]);
    % Demod and compare preamble data
    rxPreambleIndices = pamdemod(rxPreambleSymbols, 2);
    txPreambleIndices = pamdemod(txPreambleSymbols, 2);
    preErr(bn, 1, trG) = sum(txPreambleIndices ~= rxPreambleIndices);
    disp(['Preamble symbol error ' num2str(preErr(bn, 1, trG))]);

    % Demod received symbols
    rxIndices = qamdemod(rxSymbols, M);

    % Compare the payload data ignoring the ignoreN last symbols
    err(bn, :, trG) = txIndices(1:end-ignoreN) ~=...
              rxIndices(1:end-ignoreN);
    disp(['Payload symbol error ' num2str(sum(err(bn, :, trG)))]);
    
    % Compute SNR ignoring the ignoreN last symbols
    snr(bn, 1, trG) = compute_snr(...
        txSymbols(1:end-ignoreN),...
        rxSymbols(1:end-ignoreN));
    disp(['SNR ' num2str(snr(bn, 1, trG))]);

    % Plot preamble and payload symbols and eye diagram from the first 1k
    % symbols of payload
    if debug
        subplot(2, 1, 1);
        plot(rxPreambleSymbols, 'x');
        subplot(2, 1, 2);
        plot(rxSymbols, 'x');
        eyediagram(...
                   rxBlock(preambleLength*L+gDelay+1:...
                           gDelay+(preambleLength+1e3)*L),...
                   4*L,...
                   4*L);
    end
end
end
end

function snr = compute_snr(txSignal, rxSignal)
    % Normalize signal power
    rxSignal = (rms(txSignal) / rms(rxSignal)) * rxSignal;

    % Extract signal noise
    noise = txSignal - rxSignal;
    % Compute the SNR
    snr = 20 * log10(rms(txSignal) / rms(noise));
end

function  [correctedSignal, wOffset, phiOffset] = correct_offset(signal, idx)
    [wOffset, phiOffset] = fft_offset(signal(idx), 2^20, 2);
    n = 0:length(signal)-1;
    correctedSignal = signal .* exp(-1j * (wOffset*n+phiOffset));
end
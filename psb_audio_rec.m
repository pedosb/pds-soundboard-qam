[rx, rxFs] = audioread('audio_samepower.wav');
%[rx, rxFs] = audioread('audio_lowsnr.wav');
%[rx, rxFs] = audioread('audio_nobackground.wav');

rxLen = length(rx);
n = (0:(rxLen-1));

rxBb = rx' .* exp(1j .* pi * centerReplic .* n);

rxResampled = resample(rxBb, 1, wavOver);

start = 1;
while true
    blockSize = nPreambleSymbols*L+gDelay;
    block = rxResampled(start:start+blockSize-1);
    [preCorr, preCorrLags] = xcorr(block, txFilteredPreamble);
    [preCorrMax, preCorrMaxIdx] = max(abs(preCorr));
    if preCorrMax / rms(block) > 100
        break
    end
    start = start + floor(blockSize / 2);
    if start + blockSize > length(rxResampled)
        error('Could not find message singal');
    end
end

preStart = start + preCorrLags(preCorrMaxIdx);

rxAllSymbols = rxResampled(preStart+gDelay-1:L:end);
rxPreambleSymbols = rxAllSymbols(1:nPreambleSymbols);
rxPayloadSymbols = rxAllSymbols(nPreambleSymbols+1:nPreambleSymbols+nPayloadSymbols);

anDiff = mean(angle(txPreambleSymbols .* rxPreambleSymbols));
rxAdjustedPreambleSymbols = rms(pammod([0 1] , 2)) .* (rxPreambleSymbols / rms(rxPreambleSymbols)) .* exp(-1j * anDiff);
rxAdjustedPayloadsymbols = rms(qammod(0:M-1, M)) .* (rxPayloadSymbols / rms(rxPayloadSymbols)) .* exp(-1j * anDiff);

plot(rxAdjustedPreambleSymbols, 'x')
hold on
plot(rxAdjustedPayloadsymbols, 'xr')
hold off

rxPayloadIndices = qamdemod(rxAdjustedPayloadsymbols, M);

nError = sum(rxPayloadIndices ~= txPayloadIndices);
noise = txPayloadSymbols - rxAdjustedPayloadsymbols;
eSNR = 20 * log10(rms(rxAdjustedPayloadsymbols) / rms(noise));

disp(['Error ' num2str(nError) ' (' num2str(nError/length(txPayloadIndices)) '%)']);
disp(['SNR   ' num2str(eSNR)]);
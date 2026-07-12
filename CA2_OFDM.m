%% CA2_OFDM_Project.m

clear; clc; close all;
rng(1405);

%% تنظیمات اصلی پروژه
cfg.Nbits        = 1e7;
cfg.nc           = 400;
cfg.Nfft         = 1024;
cfg.M            = 4;
cfg.sfc          = ceil(2^13/cfg.nc);
cfg.cpLen        = round(0.25*cfg.Nfft);
cfg.headerFactor = 8;
cfg.snrA_dB      = 20;
cfg.snrRange_dB  = 0:2:20;
cfg.maxDelay     = [];

fprintf('Generating random bits...\n');
bits = randi([0 1], cfg.Nbits, 1, 'uint8');

fprintf('Building OFDM transmitter signal...\n');
[tx, meta] = ofdm_tx(bits, cfg);

fprintf('\nSystem parameters:\n');
fprintf('sfc = %d\n', cfg.sfc);
fprintf('CP length = %d\n', cfg.cpLen);
fprintf('OFDM symbols per frame = %d\n', cfg.sfc + 1);
fprintf('Time samples per frame = %d\n', meta.frameLen);
fprintf('Number of frames = %d\n', meta.numFrames);
fprintf('Number of padded QPSK symbols = %d\n', meta.padSymbols);
fprintf('Useful signal power = %.6g\n', meta.signalPower);

%% -------------------- بخش الف: AWGN در SNR = 20 dB --------------------
fprintf('\nPart A: AWGN, SNR = %.1f dB\n', cfg.snrA_dB);
rxA = add_awgn(tx, cfg.snrA_dB, meta.signalPower);
bitsHatA = ofdm_rx(rxA, cfg, meta);
berA = mean(bitsHatA ~= bits);

fprintf('BER Part A = %.6g\n', berA);
fprintf('Number of frames Part A = %d\n', meta.numFrames);

%% -------------------- بخش ب: نمودار BER بر حسب SNR برای AWGN --------------------
fprintf('\nPart B: BER curve for AWGN channel...\n');
berAWGN = zeros(size(cfg.snrRange_dB));

for ii = 1:numel(cfg.snrRange_dB)
    snrDb = cfg.snrRange_dB(ii);
    rx = add_awgn(tx, snrDb, meta.signalPower);
    bitsHat = ofdm_rx(rx, cfg, meta);
    berAWGN(ii) = mean(bitsHat ~= bits);
    fprintf('AWGN: SNR = %2d dB, BER = %.6g\n', snrDb, berAWGN(ii));
end

figure;
semilogy(cfg.snrRange_dB, max(berAWGN, 1/numel(bits)), '-o', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER of OFDM-QPSK-DPSK over AWGN Channel');

%% -------------------- بخش ج: نمودار BER بر حسب SNR برای کانال Rayleigh --------------------
fprintf('\nPart C: BER curve for Rayleigh channel...\n');
berRayleigh = zeros(size(cfg.snrRange_dB));

for ii = 1:numel(cfg.snrRange_dB)
    snrDb = cfg.snrRange_dB(ii);
    rx = add_rayleigh_flat_per_frame(tx, cfg, meta, snrDb);
    bitsHat = ofdm_rx(rx, cfg, meta);
    berRayleigh(ii) = mean(bitsHat ~= bits);
    fprintf('Rayleigh: SNR = %2d dB, BER = %.6g\n', snrDb, berRayleigh(ii));
end

figure;
semilogy(cfg.snrRange_dB, max(berRayleigh, 1/numel(bits)), '-s', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER of OFDM-QPSK-DPSK over Rayleigh Channel');

figure;
semilogy(cfg.snrRange_dB, max(berAWGN, 1/numel(bits)), '-o', 'LineWidth', 1.5);
hold on;
semilogy(cfg.snrRange_dB, max(berRayleigh, 1/numel(bits)), '-s', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER Comparison: AWGN vs Rayleigh');
legend('AWGN', 'Rayleigh', 'Location', 'southwest');

%% -------------------- بخش د امتیازی: همزمان سازی با Header --------------------
fprintf('\nPart D Bonus: Unknown delay and header detection...\n');

maxDelay = meta.frameLen;
[txDelayed, trueDelay] = add_random_prefix_before_header(tx, maxDelay, meta.signalPower);
rxD = add_awgn(txDelayed, cfg.snrA_dB, meta.signalPower);

searchMaxDelay = maxDelay;
detectedStart = detect_header_start(rxD, meta.header, searchMaxDelay);
rxSync = rxD(detectedStart : detectedStart + meta.totalLen - 1);

bitsHatD = ofdm_rx(rxSync, cfg, meta);
berD = mean(bitsHatD ~= bits);

fprintf('True delay = %d samples\n', trueDelay);
fprintf('Detected start index = %d\n', detectedStart);
fprintf('Detected delay = %d samples\n', detectedStart - 1);
fprintf('BER after synchronization = %.6g\n', berD);

%% ============================ Local Functions ============================

function [tx, meta] = ofdm_tx(bits, cfg)

    [symIdx, bitPad] = bits_to_qpsk_index(bits);
    numInfoSymbols = numel(symIdx);

    symbolsPerFrame = cfg.sfc * cfg.nc;
    padSymbols = mod(symbolsPerFrame - mod(numInfoSymbols, symbolsPerFrame), symbolsPerFrame);
    symIdxPadded = [symIdx; zeros(padSymbols, 1)];

    numFrames = numel(symIdxPadded) / symbolsPerFrame;
    frameLen = (cfg.Nfft + cfg.cpLen) * (cfg.sfc + 1);
    guardLen = frameLen;
    headerLen = cfg.headerFactor * frameLen;

    frames = cell(numFrames, 1);
    powerSum = 0;
    sampleCount = 0;

    for k = 1:numFrames
        st = (k-1)*symbolsPerFrame + 1;
        en = k*symbolsPerFrame;
        dataMat = reshape(symIdxPadded(st:en), cfg.sfc, cfg.nc);
        frameSig = build_one_ofdm_frame(dataMat, cfg);
        frames{k} = frameSig;
        powerSum = powerSum + sum(abs(frameSig).^2);
        sampleCount = sampleCount + numel(frameSig);
    end

    signalPower = powerSum / sampleCount;

    header = sqrt(signalPower) * (2*randi([0 1], headerLen, 1) - 1);
    guard = zeros(guardLen, 1);

    totalLen = headerLen + numFrames*(guardLen + frameLen) + guardLen + headerLen;
    tx = zeros(totalLen, 1);
    frameStarts = zeros(numFrames, 1);

    ptr = 1;
    tx(ptr:ptr+headerLen-1) = header;
    ptr = ptr + headerLen;

    for k = 1:numFrames
        tx(ptr:ptr+guardLen-1) = guard;
        ptr = ptr + guardLen;

        frameStarts(k) = ptr;
        tx(ptr:ptr+frameLen-1) = frames{k};
        ptr = ptr + frameLen;
    end

    tx(ptr:ptr+guardLen-1) = guard;
    ptr = ptr + guardLen;
    tx(ptr:ptr+headerLen-1) = header;

    meta.numOriginalBits = numel(bits);
    meta.bitPad = bitPad;
    meta.numInfoSymbols = numInfoSymbols;
    meta.padSymbols = padSymbols;
    meta.symbolsPerFrame = symbolsPerFrame;
    meta.numFrames = numFrames;
    meta.frameLen = frameLen;
    meta.guardLen = guardLen;
    meta.headerLen = headerLen;
    meta.totalLen = totalLen;
    meta.header = header;
    meta.frameStarts = frameStarts;
    meta.signalPower = signalPower;
end

function frameSig = build_one_ofdm_frame(dataMat, cfg)

    refRow = randi([0 cfg.M-1], 1, cfg.nc);
    rawIdx = [refRow; dataMat];

    dpskIdx = rawIdx;
    for r = 2:size(rawIdx, 1)
        dpskIdx(r, :) = mod(rawIdx(r, :) + dpskIdx(r-1, :), cfg.M);
    end

    qpskDpskSym = exp(1j * 2*pi/cfg.M * dpskIdx);

    D = qpskDpskSym.'; % nc x (sfc+1)
    X = zeros(cfg.Nfft, cfg.sfc + 1);
    X(2:cfg.nc+1, :) = D;
    X(cfg.Nfft-cfg.nc+1:cfg.Nfft, :) = conj(flipud(D));

    x = real(ifft(X, cfg.Nfft, 1));
    xcp = [x(end-cfg.cpLen+1:end, :); x];
    frameSig = xcp(:);
end

function bitsHat = ofdm_rx(rx, cfg, meta)

    rx = rx(:);
    symHatAll = zeros(meta.numFrames * meta.symbolsPerFrame, 1);
    outPtr = 1;
    ptr = meta.headerLen + 1;

    for k = 1:meta.numFrames
        ptr = ptr + meta.guardLen;
        frameRx = rx(ptr:ptr+meta.frameLen-1);
        ptr = ptr + meta.frameLen;

        symHat = decode_one_ofdm_frame(frameRx, cfg);
        symHatAll(outPtr:outPtr+numel(symHat)-1) = symHat;
        outPtr = outPtr + numel(symHat);
    end

    symHatAll = symHatAll(1:meta.numInfoSymbols);
    bitsHat = qpsk_index_to_bits(symHatAll);
    bitsHat = uint8(bitsHat(1:meta.numOriginalBits));
end

function symHat = decode_one_ofdm_frame(frameRx, cfg)
    ycp = reshape(frameRx, cfg.Nfft + cfg.cpLen, cfg.sfc + 1);
    y = ycp(cfg.cpLen+1:end, :);
    Y = fft(y, cfg.Nfft, 1);

    Dhat = Y(2:cfg.nc+1, :).'; % (sfc+1) x nc

    z = Dhat(2:end, :) .* conj(Dhat(1:end-1, :));

    phaseStep = 2*pi/cfg.M;
    symHatMat = mod(round(angle(z) / phaseStep), cfg.M);
    symHat = symHatMat(:);
end

function [symIdx, bitPad] = bits_to_qpsk_index(bits)
    bits = uint8(bits(:));
    bitPad = 0;
    if mod(numel(bits), 2) ~= 0
        bits = [bits; uint8(0)];
        bitPad = 1;
    end

    symIdx = double(bits(1:2:end))*2 + double(bits(2:2:end));
end

function bits = qpsk_index_to_bits(symIdx)
    symIdx = double(symIdx(:));
    bits = zeros(2*numel(symIdx), 1, 'uint8');
    bits(1:2:end) = uint8(floor(symIdx/2));
    bits(2:2:end) = uint8(mod(symIdx, 2));
end

function rx = add_awgn(tx, snrDb, signalPower)
    snrLin = 10^(snrDb/10);
    noiseVar = signalPower / snrLin;
    rx = tx + sqrt(noiseVar) * randn(size(tx));
end

function rx = add_rayleigh_flat_per_frame(tx, cfg, meta, snrDb)
    rx = tx;
    U = max(rand(meta.numFrames, 1), realmin);
    h = sqrt(-log(U));

    for k = 1:meta.numFrames
        idx = meta.frameStarts(k) : meta.frameStarts(k) + meta.frameLen - 1;
        rx(idx) = h(k) * rx(idx);
    end

    rx = add_awgn(rx, snrDb, meta.signalPower);
end

function [txDelayed, delay] = add_random_prefix_before_header(tx, maxDelay, signalPower)
    delay = randi([0 maxDelay], 1, 1);
    prefix = sqrt(signalPower) * randn(delay, 1);
    txDelayed = [prefix; tx(:)];
end

function startIdx = detect_header_start(rx, header, searchMaxDelay)

    rx = rx(:);
    header = header(:);
    M = numel(header);

    searchLen = min(numel(rx), searchMaxDelay + M);
    x = rx(1:searchLen);
    N = numel(x);

    L = 2^nextpow2(N + M - 1);
    c = ifft(fft(x, L) .* fft(flipud(conj(header)), L));
    r = c(M:N);

    cs = [0; cumsum(abs(x).^2)];
    winEnergy = cs(M+1:end) - cs(1:end-M);
    headerEnergy = sum(abs(header).^2);

    metric = abs(r).^2 ./ (winEnergy * headerEnergy + eps);
    [~, startIdx] = max(metric);
end

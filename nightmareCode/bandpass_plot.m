function bandpass_plot(signal, Fs)
    % BANDPASS_PLOT filters a signal between 0.5 and 3 Hz,
    % plots the power spectrum (DFT), and plots the spectrogram.
    %
    % signal : input vector
    % Fs     : sampling frequency (Hz)

    % --- Design bandpass filter ---
    bpFilt = designfilt('bandpassiir', ...
                        'FilterOrder', 4, ...
                        'HalfPowerFrequency1', 0.5, ...
                        'HalfPowerFrequency2', 3, ...
                        'SampleRate', Fs);

    % --- Apply zero-phase filter ---
    filtered_signal = filtfilt(bpFilt, signal);

    % --- Compute DFT ---
    N = length(filtered_signal);
    X = fft(filtered_signal);
    f = (0:N-1)*(Fs/N);
    powerX = abs(X).^2 / N;

    figure;
    subplot(3,1,1)
    plot((1:length(signal))./Fs, signal)
    xlabel('Time (s)')
    ylabel('ROI Pixel Sum')
    xlim([1 15])

    % --- Plot Power Spectrum ---
    
    subplot(3,1,2)
    plot(f, powerX);
    xlim([0 10]);  % Only low-frequency content matters
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title('Power Spectrum of Band-Passed Signal');
    grid on;

    % --- Plot Spectrogram ---
    subplot(3,1,3)
    window = hamming(round(2*Fs)); % 2-second window
    overlap = round(0.9 * length(window)); % 90% overlap
    nfft = 1024;

    spectrogram(filtered_signal, window, overlap, nfft, Fs, 'yaxis');
    colormap jet;
    title('Spectrogram of Band-Passed Signal');
    ylim([0 5]);   % Only show low frequencies of interest
end
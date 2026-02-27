function respirationSpectrogramPlot(signal, Fs)
    % BANDPASS_PLOT filters a signal between 0.5 and 3 Hz,
    % plots the power spectrum (DFT), and plots the spectrogram
    % with frequency center of mass overlay.
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
    subplot(4,1,1)
    plot((1:length(signal))./Fs, signal)
    xlabel('Time (s)')
    ylabel('ROI Pixel Sum')
    title('Respiration ROI (Front Chest) Pixel Sum');
    xlim([1 15])

    subplot(4,1,2)
    plot((1:length(filtered_signal))./Fs, filtered_signal)
    xlabel('Time (s)')
    ylabel('Filtered ROI Pixel Sum')
    title('Band-Passed (0.5-3Hz) Respiration ROI (Front Chest) Pixel Sum');
    xlim([1 15])

    % --- Plot Power Spectrum ---
    subplot(4,1,3)
    plot(f, powerX);
    xlim([0 10]);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title('Power Spectrum of Band-Passed Signal');
    grid on;

    % --- Plot Spectrogram ---
    subplot(4,1,4)

    window  = hamming(round(2*Fs));      % 2-second window
    overlap = round(0.9 * length(window));
    nfft    = 1024;

    % Compute spectrogram explicitly (so we can use the data)
    [S,F,T] = spectrogram(filtered_signal, window, overlap, nfft, Fs);

    % Power spectrogram
    P = abs(S).^2;

    % Plot spectrogram
    imagesc(T, F, 10*log10(P));
    axis xy
    colormap jet;
    ax = gca;
    cb = colorbar(ax, 'eastoutside');
    ax.PositionConstraint = 'innerposition';
    ylabel(cb, 'Power')
    title('Spectrogram of Band-Passed Signal');
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    ylim([0 5])
    xlim([1 15])

    % --- CENTER OF MASS CALCULATION (Frequency Centroid) ---
    % Weighted mean frequency at each time slice
    freq_centroid = sum(F .* P, 1) ./ sum(P, 1);

    % --- Overlay center of mass ---
    hold on;
    plot(T, freq_centroid, 'w', 'LineWidth', 2);
    plot(T, freq_centroid, 'k--', 'LineWidth', 1); % outline for contrast
    hold off;
end

function [tvec, Evec] = tACSwave(dt, DEL, DUR, tstop, waveform, plot_wave)
    amp = waveform.amp;
    freq = waveform.freq;

    tvec_del = (0:dt:DEL-dt)'; 
    tvec_dur = (0:dt:DUR-dt)'; 
    tvec_tstop = (0:dt:tstop-DEL-DUR)';
    tvec = [tvec_del; tvec_dur + DEL; tvec_tstop + DEL + DUR];
    
    Evec_del = zeros(size(tvec_del), 'like', tvec_del);
    Evec_dur = amp * sin(2*pi*freq*tvec_dur/1000);
    Evec_tstop = zeros(size(tvec_tstop), 'like', tvec_tstop);
    Evec = [Evec_del; Evec_dur; Evec_tstop];

    if plot_wave
        figure('Color','w');
        hold on;
        plot(tvec,Evec,'Color','k','LineWidth',2)
        title('tACS waveform');
        xlabel('time (ms)'); xlim([0, tstop]);
        ylabel('Amplitude of tACS (mV/mm)');
        box off;
        set(gca,'FontSize',20)
    end
end
function closeVDAM(measurement_angle, totalData,filepath,MA)
    fprintf('saving data, closing VDAM...\n');
    filename = [filepath '\' datestr(now,'yyyy-mm-dd_HHMMSS') '_' measurement_angle, '.mat'];
    save(filename, 'totalData');
    MA.stop();
    MA.close();
    delete(MA);
end
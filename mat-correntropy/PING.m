function z = PING (host, times)
    cmd = ['!ping -c ',times,' ', host];
    strping = evalc(cmd);
    z = cellfun(@(y)(str2double(y{2})), regexp(strping, '(time[=<])([\d]*)', 'tokens'));
end
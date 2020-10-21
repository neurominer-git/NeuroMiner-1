function str = timestampstr
str = char(datetime(datestr(now),'Format','ddMMyyyy_HHmmss'));
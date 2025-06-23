%此函数用来将计算出来的分和秒用字符串的形式进行表示
function minute_and_second = get_minute_and_second (k)
minute = addLeadingZeroIfSingleDigit(floor(k / 2));
second = addLeadingZeroIfSingleDigit(30*mod(k,2));
minute_and_second = strcat(minute, ':', second);

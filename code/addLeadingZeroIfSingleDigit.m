%此函数用来在单数前面加一个0变成字符串
function formattedNumber = addLeadingZeroIfSingleDigit(number)
    % 将数字转换为字符串
    numStr = num2str(number);
    % 检查数字是否是个位数
    if length(numStr) == 1 && numStr ~= '-'
        % 如果是，添加一个0到字符串前面
        formattedNumber = ['0', numStr];
    else
        % 如果不是，直接返回原字符串
        formattedNumber = numStr;
    end
end
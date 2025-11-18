% black box manual method makeshift func


blackBoxFunc = @griewank3;

disp("Welcome to the griewank function");

while true
    point = [];
    x1 = loopTillNum("Enter x1: ");
    x2 = loopTillNum("Enter x2: ");
    x3 = loopTillNum("Enter x3: ");

    point = [x1, x2, x3];
    result = blackBoxFunc(point);
    disp(['Function output: ', num2str(result)]);



end

function x = loopTillNum(inputMessage)

    x = nan;
    while isnan(x)
            userInput = input(inputMessage, 's');
            x = str2double(userInput);
            if isnan(x)
                disp('Invalid input. Please enter a numeric value.');
            end
    end


end
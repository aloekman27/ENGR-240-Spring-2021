yardLine = 10;
v0 = 65;
theta0 = 35;
[missedByOrMadeBy, points] = fieldGoal(yardLine, v0, theta0)

function [missedByOrMadeBy, points] = fieldGoal(yardLine, v0, theta0)
%generates a constant where y0 = 0ft and g = 32.2ft/s^2 from the intial problem statement
y0 = 0;
g = 32.2;
%an if statement for when the number of elements or the type of vector of the variable yardline & v0 does not match. It will produce an output of 999.
if length(yardLine) ~= length(v0)
    missedByOrMadeBy= 999;
    points = 999;
%create a column vector of zeros from the number of elements in yardline. This will be useful for the output in order to be a column vector
else missedByOrMadeBy = zeros(length(yardLine),1);
    %each point that is gotten will be replacing the corresponding zeros in the column vector
    points = missedByOrMadeBy;
    %use a for loop to determine how many times MATLAB should redo the calculation steps
    for s = 1:length(v0) %number of loops is up to the number of elements in the variable yardLine (or v0)
        x = (yardLine(s) + 17) * 3; % each x is equal to the yardLine (in ft) added with 17 (from the end of the football field) and multiplied by 3 (convert yards into feet)
        y = x*tand(theta0) - ((1/2)*(x^2*g)./((v0(s)*cosd(theta0))^2)) + y0; % calculate the final height using trajectory formula
        missedByOrMadeBy(s) = y-10; %the final height will be subtracted by 10 in order to know whether the kicker makes the field goal or not. Also, this will be shown in the output
        
        %use if statements to determine the number of points gotten for different situations
        if missedByOrMadeBy(s) < 0
            points(s) = 0; %if the final height subtracted by 10 has a negative value, it will receive 0 points since it means that the kicker does not make the field goal
        else if yardLine(s) == 15
                points(s) = 1; % special case wehre the yardLine = 15. This means that the team is attempting for an extra point if successful (final height - 10 > 0)
        else points(s) = 3; % all other cases except for the mentioned above will receive 3 points
        end 
        end 
    end 
end
end
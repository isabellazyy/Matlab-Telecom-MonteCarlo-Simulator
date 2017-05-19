function ul = user_location(numOfUsers,numOfIterations,meter)

    ul = zeros(numOfUsers * numOfIterations,2);

    for i = 1 : numOfIterations * numOfUsers
        ul(i,:) = randi([0 meter],1,2);
    end
    
end


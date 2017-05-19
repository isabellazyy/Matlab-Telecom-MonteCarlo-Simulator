function d = distance(numOfUsers,numOfCells,numOfIterations,meter)

    C = cell_location(numOfCells,meter);
    U = user_location(numOfUsers,numOfIterations,meter);
    d = zeros(numOfCells, numOfIterations *numOfUsers);
    
    for i = 1:numOfCells
        for k = 1:numOfUsers * numOfIterations
            temp = C(i,:) - U(k,:);
            % distance between users and cells
            d(i,k) = sqrt(sum((temp.^2)));   
        end 
    end
end
function d = distance(numOfUsers,numOfCells,numOfIterations,meter)

    CL = randi([0 meter],numOfCells,2);
 
    d = zeros(numOfCells,numOfUsers,numOfIterations);
     
    for i=1:numOfIterations
   
       U = randi([0 meter],numOfUsers,2);
        
       for j = 1:numOfCells
                  
           for k = 1:numOfUsers
            temp = CL(j,:) - U(k,:);
            % distance between users and cells
            d(j,k,i) = sqrt(sum((temp.^2))); 
               
           end
       end
           
    end
    
end
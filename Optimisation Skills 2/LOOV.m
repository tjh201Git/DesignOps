function RMSE_LOO = LOOV(Nodes,Phi,MCSI_samples)
yhatLOO = zeros(Nodes,1); %Initialize leave one out yhat
sseLOO = 0;              %Initialize sum of squared errors  

for i = 1:Nodes 

    %CROSS VALIDATION - use the data of all nodes except the current one
    idx_train = setdiff(1:Nodes, i);  
    
    Phi_train = Phi(idx_train,:);  %trains a phi using the current points except the i'th point
    y_train   = MCSI_samples(idx_train);      
    A_i = pinv(Phi_train) * y_train;
    yhat_i = Phi(i,:) * A_i;   %Makes a surrogate model using every point except i
                               %a different surrogate is made each time in
                               %loop
    yhatLOO(i) = yhat_i;       %Stores the current yhat leave one out value
    
    
    sseLOO = sseLOO + (MCSI_samples(i) - yhat_i)^2;  %calculates the total sum of squared error
                                                     % Accumulate error
                                                     % across all
end

RMSE_LOO = sqrt(sseLOO / Nodes); %Root mean squared of the sum of squared errors

end

for j=2:loop            % Time-loop
    {
    for i=1:len(r)-1    %Space-loop
        {
            t=j*dt;
            u=(1/4).*(1./(1-t/tf)).*1./(r(i)/R).*((1-(r(i)/R.^2)).^-lambda - (1-(r(i)/R).^2));
         }
    }
    end
    
end

for j=2:loop 
           
           t=j*dt;
           u=(1/4).*(1./(1-t/tf)).*1./(r(1:end-1)./R).*((1-(r(1:end-1)./R.^2)).^-lambda - (1-(r(1:end-1)./R).^2));
end
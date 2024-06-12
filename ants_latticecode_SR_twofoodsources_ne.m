clear

dt = .001; % Time Step
dx =  1.0;  % Spatial Step
N = 50;     % Lattice Size N+1 x N+1
M = 1000;     % Number of Ants
Mf = 2;     % Number of Food Sources
D = 10.0;     % Diffusion
gamma = .001;    % Evaporation Coefficient
CFL = D*dt/(dx*dx)  % CFL Condition
q = 1.0;         % Ant Pheromone Production Rate
Time = 2.0;
count1 = 0;
Track = 0;

% Build Lattice
for i = 1:N+1
    for j = 1:N+1
        x(i,j) = (i-1)*dx;
        y(i,j) = (j-1)*dx;
    end
end

%IC on Pheromone
for i = 1:N+1
    for j = 1:N+1
        c(i,j,1) = 0.0;
        c(i,j,2) = 0.0;
        dens(i,j) = 0.0;
        dens_ret(i,j) = 0.0;
    end
end

%Place Ants at Nest y
xn(1) = 25;
xn(2) = 25;

for i = 1:M
    ant(i,1) = xn(1);
    ant(i,2) = xn(2); %Place nest at center
end
% Food Source Location
    xf(1,1) = 5.0;
    xf(1,2) = 45.0;
    xf(2,1) = 25.0;
    xf(2,2) = 5.0;

food_found = 0; %changes to 1 when food is found
for i = 1:M
    mark(i) = 0; %0 for foragers, 1 returners
end

% Time Loop
for k = 1:Time/dt
        % Check if Food is Found
        for p = 1:M
            for p1 = 1:Mf
            if (ant(p,1) == xf(p1,1) && ant(p,2) == xf(p1,2))
                food_found = 1;
            end
            end
        end
        % Check conversion for forager to returner
        for i = 1:M
            if (mark(i) == 0)
                for p1 = 1:Mf
                    dist1 = sqrt((ant(i,1)-xf(p1,1))^2+(ant(i,2)-xf(p1,2))^2);
                      if (dist1 <= dx)
                      mark(i) = 1;  % change to returner
                      fs(i,1) = xf(p1,1);
                      fs(i,2) = xf(p1,2);
                      end
                end
            elseif (mark(i) == 1)
                 dist1 = sqrt((ant(i,1)-xn(1))^2+(ant(i,2)-xn(2))^2);
                      if (dist1 <= dx)
                      mark(i) = 0;  % change to forager
                      end
            end
        end
        
        if (food_found > 0)
            %Solve PDE for c(x,t) with BCs
 %       'Food Found'
            for i = 2:N
                for j = 2:N
                    c(i,j,2) = c(i,j,1)+(D*dt/(dx*dx))*(c(i+1,j,1)+c(i-1,j,1)-4.0*c(i,j,1)+c(i,j+1,1)+c(i,j-1,1))-gamma*dt*c(i,j,1);
                end
            end
            % Add Ant Pheromone Contribution
            for p = 1:M
                if (mark(p) == 1)
                    c(ant(p,1),ant(p,2),2) = c(ant(p,1),ant(p,2),2) + dt*q*exp(-((ant(p,1)-fs(p,1))^2+(ant(p,2)-fs(p,2))^2));
                end
            end
            % BCs
            for i = 1:N+1
                c(1,i,2) = c(2,i,1)-dx/D*c(1,i,1);
                c(N+1,i,2) = c(N,i,1)-dx/D*c(N+1,i,1);
                c(i,1,2) = c(i,2,1)-dx/D*c(i,1,1);
                c(i,N+1,2) = c(i,N,1)-dx/D*c(i,N+1,1);
            end

            % Reset for next time step
                for i = 1:N+1
                    for j = 1:N+1
                        if (c(i,j,2) < 0.000000000001)
                            c(i,j,2) = 0; % Ensure pheromone concnetration is nonnegative
                        end
                        c(i,j,1) = c(i,j,2);
                    end
                end
        end


    for i = 1:M
        if (mark(i) == 1) % returner dynamics
        vec_home(1) = xn(1)-ant(i,1);
        vec_home(2) = xn(2)-ant(i,2);
        
        tempx = ant(i,1) + 1.0*vec_home(1)/sqrt(vec_home(1)^2+vec_home(2)^2);
        tempy = ant(i,2) + 1.0*vec_home(2)/sqrt(vec_home(1)^2+vec_home(2)^2);

        % % Find the nearest Grid Point
        % indx = tempx/dx+1;
        % indy = tempy/dx+1;
        % dist1 = 10.0;
        % % New Point Candidate 1
        % dist2 = sqrt((tempx-indx)^2+(tempy-indy)^2);
        % if (dist2 < dist1)
        %     dist1 = dist2;
        %     i1 = indx;
        %     i2 = indy;
        % end
        % % New Point Candidate 2
        % if (vec_home(1) < 0)
        % dist2 = sqrt((tempx-(indx-1))^2+(tempy-indy)^2);
        % if (dist2 < dist1)
        %     dist1 = dist2;
        %     i1 = indx-1;
        %     i2 = indy;
        % end
        % elseif (vec_home(1) > 0)
        % dist2 = sqrt((tempx-(indx+1))^2+(tempy-indy)^2);
        % if (dist2 < dist1)
        %     dist1 = dist2;
        %     i1 = indx+1;
        %     i2 = indy;
        % end
        % end
        % % New Point Candidate 3
        % if (vec_home(2) < 0)
        % dist2 = sqrt((tempx-indx)^2+(tempy-(indy-1))^2);
        % if (dist2 < dist1)
        %     dist1 = dist2;
        %     i1 = indx;
        %     i2 = indy-1;
        % end
        % elseif (vec_home(2) > 0)
        % dist2 = sqrt((tempx-indx)^2+(tempy-(indy+1))^2);
        % if (dist2 < dist1)
        %     dist1 = dist2;
        %     i1 = indx;
        %     i2 = indy+1;
        % end
        % end
        % % New Point Candidate 4
        % if (vec_home(1) < 0 && vec_home(2) < 0)
        % dist2 = sqrt((tempx-(indx-1))^2+(tempy-(indy-1))^2);
        % if (dist2 < dist1)
        %     dist1 = dist2;
        %     i1 = indx-1;
        %     i2 = indy-1;
        % end
        % elseif (vec_home(1) > 0 && vec_home(2) > 0)
        % dist2 = sqrt((tempx-(indx+1))^2+(tempy-(indy+1))^2);
        % if (dist2 < dist1)
        %     dist1 = dist2;
        %     i1 = indx+1;
        %     i2 = indy+1;
        % end
        % end


        ant(i,1) = int16(tempx);
        ant(i,2) = int16(tempy);
       
        end
        %forager dynamics
        if (mark(i) == 0) 
        % Compute probability for 8 neighboring cells
        %Neighbor 1
        local_c = c(ant(i,1)-1,ant(i,2)+1,2)-c(ant(i,1),ant(i,2),2);
        if (local_c < 0)
            wa(i,1) = 0;
        elseif (local_c == 0)
            wa(i,1) = 1;
        else
            wa(i,1) = 1+local_c;
        end
        %Neighbor 2
        local_c = c(ant(i,1),ant(i,2)+1,2)-c(ant(i,1),ant(i,2),2);
        if (local_c < 0)
            wa(i,2) = 0;
        elseif (local_c == 0)
            wa(i,2) = 1;
        else
            wa(i,2) = 1+local_c;
        end
                %Neighbor 3
        local_c = c(ant(i,1)+1,ant(i,2)+1,2)-c(ant(i,1),ant(i,2),2);
        if (local_c < 0)
            wa(i,3) = 0;
        elseif (local_c == 0)
            wa(i,3) = 1;
        else
            wa(i,3) = 1+local_c;
        end
                %Neighbor 4
        local_c = c(ant(i,1)+1,ant(i,2),2)-c(ant(i,1),ant(i,2),2);
        if (local_c < 0)
            wa(i,4) = 0;
        elseif (local_c == 0)
            wa(i,4) = 1;
        else
            wa(i,4) = 1+local_c;
        end
                %Neighbor 5
        local_c = c(ant(i,1)+1,ant(i,2)-1,2)-c(ant(i,1),ant(i,2),2);
        if (local_c < 0)
            wa(i,5) = 0;
        elseif (local_c == 0)
            wa(i,5) = 1;
        else
            wa(i,5) = 1+local_c;
        end
                %Neighbor 6
        local_c = c(ant(i,1),ant(i,2)-1,2)-c(ant(i,1),ant(i,2),2);
        if (local_c < 0)
            wa(i,6) = 0;
        elseif (local_c == 0)
            wa(i,6) = 1;
        else
            wa(i,6) = 1+local_c;
        end
                %Neighbor 7
        local_c = c(ant(i,1)-1,ant(i,2)-1,2)-c(ant(i,1),ant(i,2),2);
        if (local_c < 0)
            wa(i,7) = 0;
        elseif (local_c == 0)
            wa(i,7) = 1;
        else
            wa(i,7) = 1+local_c;
        end
                %Neighbor 8
        local_c = c(ant(i,1)-1,ant(i,2),2)-c(ant(i,1),ant(i,2),2);
        if (local_c < 0)
            wa(i,8) = 0;
        elseif (local_c == 0)
            wa(i,8) = 1;
        else
            wa(i,8) = 1+local_c;
        end
            wa_sum(i) = sum(wa(i,:));
        %Choose random number from 0 to 1
      %  if (wa_sum(i) ~= 0)
        r1 = rand();
        if (r1 < wa(i,1)/wa_sum(i))
            % Move to neighbor 1
            ant(i,1) = ant(i,1)-1;
            ant(i,2) = ant(i,2)+1;
        elseif (r1 < (wa(i,1)+wa(i,2))/wa_sum(i))
            %move to neighbor 2
            ant(i,1) = ant(i,1);
            ant(i,2) = ant(i,2)+1;
        elseif (r1 < (sum(wa(i,1:3)))/wa_sum(i))
            %move to neighbor 3
            ant(i,1) = ant(i,1)+1;
            ant(i,2) = ant(i,2)+1;
        elseif (r1 < (sum(wa(i,1:4)))/wa_sum(i))
            %move to neighbor 4
            ant(i,1) = ant(i,1)+1;
            ant(i,2) = ant(i,2);
        elseif (r1 < (sum(wa(i,1:5)))/wa_sum(i))
            %move to neighbor 5
            ant(i,1) = ant(i,1)+1;
            ant(i,2) = ant(i,2)-1;
        elseif (r1 < (sum(wa(i,1:6)))/wa_sum(i))
            %move to neighbor 6
            ant(i,1) = ant(i,1);
            ant(i,2) = ant(i,2)-1;
        elseif (r1 < (sum(wa(i,1:7)))/wa_sum(i))
            %move to neighbor 7
            ant(i,1) = ant(i,1)-1;
            ant(i,2) = ant(i,2)-1;
        elseif (r1 <= (sum(wa(i,1:8)))/wa_sum(i))
            %move to neighbor 8
            ant(i,1) = ant(i,1)-1;
            ant(i,2) = ant(i,2);
        end
        % elseif (wa_sum == 0)
        % % Choose Random Direction of Motion
        % r1 = rand();
        % if (r1 < .125)
        %     % Move to neighbor 1
        %     ant(i,1) = ant(i,1)-1;
        %     ant(i,2) = ant(i,2)+1;
        % elseif (r1 < .250)
        %     %move to neighbor 2
        %     ant(i,1) = ant(i,1);
        %     ant(i,2) = ant(i,2)+1;
        % elseif (r1 < .375)
        %     %move to neighbor 3
        %     ant(i,1) = ant(i,1)+1;
        %     ant(i,2) = ant(i,2)+1;
        % elseif (r1 < .500)
        %     %move to neighbor 4
        %     ant(i,1) = ant(i,1)+1;
        %     ant(i,2) = ant(i,2);
        % elseif (r1 < .625)
        %     %move to neighbor 5
        %     ant(i,1) = ant(i,1)+1;
        %     ant(i,2) = ant(i,2)-1;
        % elseif (r1 < .750)
        %     %move to neighbor 6
        %     ant(i,1) = ant(i,1);
        %     ant(i,2) = ant(i,2)-1;
        % elseif (r1 < .875)
        %     %move to neighbor 7
        %     ant(i,1) = ant(i,1)-1;
        %     ant(i,2) = ant(i,2)-1;
        % elseif (r1 <= 1.000)
        %     %move to neighbor 8
        %     ant(i,1) = ant(i,1)-1;
        %     ant(i,2) = ant(i,2);
        % end
        % 
        % end
        % Reflecting BCs
        if (ant(i,1) == 1) 
            ant(i,1) = 2;
        end
        if (ant(i,1) == N+1)
            ant(i,1) = N;
        end
        if (ant(i,2) == 1) 
            ant(i,2) = 2;
        end
        if (ant(i,2) == N+1)
            ant(i,2) = N;
        end
        end
    end

  % Plot every so many frame
  count1 = count1 + 1;
  if (count1 > 9)
      count1 = 0;
    % Density plot
    for i = 1:N+1
        for j = 1:N+1
            dens(i,j) = 0.0;
            dens_ret(i,j) = 0.0;
        end
    end

       for p = 1:M
           if (mark(p) == 0)
           dens(ant(p,1),ant(p,2)) = dens(ant(p,1),ant(p,2))+1;
           elseif (mark(p) == 1)
           dens_ret(ant(p,1),ant(p,2)) = dens_ret(ant(p,1),ant(p,2))+1;
           end
       end

       %Scaling for Images
    for i = 1:N+1
        for j = 1:N+1
            % if (dens(i,j) > 3)
            % dens(i,j) = 3.0;
            % end
            % if (dens_ret(i,j) > 10)
            % dens_ret(i,j) = 10.0;
            % end
        % Plot Pheromone Concentration
                c1(i,j) = c(i,j,2);

        end
    end


      % hold on
     %  plot(xf(1,1),xf(1,2),'go')
    subplot(1,3,1);
  
       h1 = heatmap(dens')
       h1.GridVisible = 'off';
       h1.YDisplayData = flipud(h1.YDisplayData); 
      XLabels = 1:N+1;
% Convert each number in the array into a string
CustomXLabels = string(XLabels);
% Replace all but the fifth elements by spaces
CustomXLabels(mod(XLabels,5) ~= 0) = " ";

      h1.XDisplayLabels = CustomXLabels;
      h1.YDisplayLabels = flip(CustomXLabels);
       title('Forager Distribution')
       ax1 = gca
       colorbar
       colormap(ax1,hot)


    subplot(1,3,2)
       h2 = heatmap(dens_ret')
       h2.GridVisible = 'off';
       h2.YDisplayData = flipud(h2.YDisplayData); 
      XLabels = 1:N+1;
% Convert each number in the array into a string
CustomXLabels = string(XLabels);
% Replace all but the fifth elements by spaces
CustomXLabels(mod(XLabels,5) ~= 0) = " ";

      h2.XDisplayLabels = CustomXLabels;
      h2.YDisplayLabels = flip(CustomXLabels);
       title('Returner Distribution')
       colorbar
       colormap(sky)

 subplot(1,3,3)
       h = heatmap(c1(:,:)')
       h.GridVisible = 'off';
       h.YDisplayData = flipud(h.YDisplayData); 
      XLabels = 1:N+1;
% Convert each number in the array into a string
CustomXLabels = string(XLabels);
% Replace all but the fifth elements by spaces
CustomXLabels(mod(XLabels,5) ~= 0) = " ";

      h.XDisplayLabels = CustomXLabels;
      h.YDisplayLabels = flip(CustomXLabels);
    title('Pheromone Distribution')
       colorbar
       colormap(sky)      
    drawnow

        Track = Track + 1;
        filename=sprintf('forager_dens_%5.5d.dat',Track);
        dlmwrite(filename,dens,'delimiter','\t','precision',5)
        filename=sprintf('returner_dens_%5.5d.dat',Track);
        dlmwrite(filename,dens_ret,'delimiter','\t','precision',5)
        filename=sprintf('pheromone_dens_%5.5d.dat',Track);
        dlmwrite(filename,c1,'delimiter','\t','precision',5)        
      %  out=sprintf('chrom_mov_%5.5d.jpg',Track);
      %  filename=sprintf('chrom_mov_%5.5d.jpg',Track);
      %  saveas(gcf,filename)
  end

end



%function [xmin,ymin,f] = SA


    
    sigma = 1;
    temp_init = 10000;
    temp = temp_init;
    xy(1,1) = 2;
    xy(1,2) = 2;
    na = 0;
    
    while temp > 1
        tvec = abs (temp_init - temp) + 1; % reflect temp do count up
        
        % generate random states
        proposal = [normrnd(xy(tvec, 1), sigma), normrnd(xy(tvec, 2), sigma)];
        % Calculate delta abs
        deltaE = abs (g(xy(tvec,1), xy(tvec,2)) - g(proposal(1), proposal(2)));
        
        if deltaE > 0
            accpt = true;
        else 
            % calculate cost
            cost_arg = exp( ( g(xy(tvec,1), xy(tvec,2)) - g(proposal(1), proposal(2)) ) / temp );
            
            if (rand > cost_arg) % accept state
                accpt = true;
            else
                accpt = false;
            end
        end
        if accpt == true
            xy(tvec+1,1) = proposal(1);
            xy(tvec+1,2) = proposal(2);
            na = na + 1;
        end
        
        temp = temp - 1;
    end
    
    xmin = xy(tvec,1);
    ymin = xy(tvec,2);
    f = g(xy(tvec,1), xy(tvec,2));


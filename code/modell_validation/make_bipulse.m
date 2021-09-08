function [t,sig] = make_bipulse(amp, duration, t_delay, t_up, t_pause, t_down, T, ...
    delta, offset)
    
    t = delta:delta:duration;
    sig = -offset*ones(1,length(t));
    
    n = 1;
    for j = 1:length(t)
        if j > n*T/delta
            n=n+1;
        end
        lim = (n-1)*T/delta + t_delay/delta;
        if j > lim && j <= lim + t_up/delta
            sig(j) = amp;
        elseif j > lim + (t_up+t_pause)/delta && j <= lim + (t_up+t_pause+t_down)/delta
            sig(j) = -amp;
        end
    end
    sig = sig + offset;

end
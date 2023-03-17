function T = get_temperature(t, target_temperature, init_temperature, varargin)



if ~isempty(varargin)
    i = varargin{1}(1);
    j = varargin{1}(2);
    i_max = varargin{2}(1);
    j_max = varargin{2}(2);
else
    i = 1;
    j = 1;
    i_max = 0;
    j_max = 0;
end



if i > i_max && j > j_max

    if t < 5000
        T = init_temperature;
    elseif t < 50
        T = init_temperature +  (target_temperature - init_temperature)/20 * (t - 30);
%         T = target_temperature;
    elseif t < 300
        T = target_temperature;
    elseif t < 320
        T = target_temperature - (target_temperature - init_temperature)/20 * (t - 300);
    else
        T = init_temperature;
    end


else
    T = 15;
end

end
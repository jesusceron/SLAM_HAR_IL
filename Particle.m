classdef Particle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        W
        X
        Y
        Lm
        LmP
        T_x
        T_y
    end
    
    methods
        function obj = Particle(N_PARTICLES,N_LM,LM_SIZE,x,y)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.W = 1/N_PARTICLES;
            obj.X = x;
            obj.Y = y;
            obj.Lm = zeros(N_LM, LM_SIZE);
            obj.LmP = {};
            for i=1:N_LM
                obj.LmP{i} = zeros(LM_SIZE, LM_SIZE);
            end
            obj.T_x = x;
            obj.T_y = y;            
            
        end

    end
end


classdef IVP < matlab.mixin.SetGet 
    %A simple class to represent an initial value problem
    % - main ouput object from XPP ode file parser
    %
    % - provide methods for generating grids or samplings of
    % parameter/initial condition spaces 
    % - provide a simple GUI to modify parameters and initial conditions
    
    %todo: make it MatCont compatible too
    properties
    
        mfun %filename of matlab function that computes the vector field
        
        pnames %allow changes to names for graphing
        p0 %parameter set
        
        xnames %allow changes to names for graphing
        x0 %initial condition
        
        tscale=1
        tunit='' %for display
        
    end
    
    properties ( Access = private )
        
        
    end
    
    methods
        
        function ivp=IVP(varargin)
            %none - uigetfile to find ODE file
            %char - ODE filename
            %struct - parseODEfile output
        end
        
    end
    
    methods ( Access = private )
    end

end
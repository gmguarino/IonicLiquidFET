classdef Constants
    properties( Constant = true )
        %ionic liquid dielectric constant
        K_IL = 15;
        %nanocapacitor gap size
        d_nc = 1e-8;
        %width of channel
        W = 10e-6;
        %length of the channel
        L= 100e-6;
        
        
        %elementary charge
        e = - 1.602e-19;
        epsilon_0 = 8.854e-12;
        hbar = 1.05457182e-34;
        m_e = 9.11e-31;
        %dielectric constant of STO
        K_STO = 300;
        
        ml = 1.2 .* m_e;
        gl = 4;
        mh = 4.8 .* m_e;
        gh = 2;

    end
 end
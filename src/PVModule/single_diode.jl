"""
Base model definition for the single diode model.

For a detailed description of the model, see:
https://pvpmc.sandia.gov/modeling-guide/2-dc-module-iv/single-diode-equivalent-circuit-models/

The diode light current is specified using callable parameters:
https://docs.sciml.ai/ModelingToolkit/dev/tutorials/callable_params/#Callable-parameters-and-interpolating-data
"""

@mtkmodel HeatingDiode begin
    begin
        k = 1.380649e-23 # Boltzmann constant (J/K)
        q = 1.602176634e-19 # Elementary charge (C)
    end

    @extend v, i = oneport = OnePort(; v = 0.0)
    @components begin
        port = HeatPort()
    end
    @parameters begin
        Is = 1e-6, [description = "Saturation current (A)"]
        n = 1, [description = "Ideality factor"]
    end
    @variables begin
        Vt(t), [description = "Thermal voltage"]
    end
    @equations begin
        Vt ~ k * port.T / q  # Thermal voltage equation
        i ~ Is * (exp(v / (n * Vt)) - 1)  # Shockley diode equation
        port.Q_flow ~ -v * i  # -LossPower
    end
end

function SingleDiode(data, time)
    @named I_L = ParametrizedInterpolation(LinearInterpolation, data, time)

    @mtkmodel Single_Diode begin
        @parameters begin
            n::Float64 = 1.0, [description = "Ideality factor"]
            I_s::Float64 = 1e-6, [description = "Saturation current (A)"]
            R_sh_val::Float64 = 1e5, [description = "Shunt resistance (Ohm)"]
            R_s_val::Float64 = 1e-1, [description = "Series resistance (Ohm)"]
        end

        @components begin
            R_sh = Resistor(R = R_sh_val)
            R_s = Resistor(R = R_s_val)
            diode = HeatingDiode(n = n, Is = I_s)
            source = Current()
            temp = FixedTemperature(T = 300.0)
            pos = Pin()
            neg = Pin()
        end

        @equations begin

            # current from irradiance
            connect(I_L.output, source.I)

            # positive side
            connect(source.n, diode.p)
            connect(source.n, R_sh.p)
            connect(source.n, R_s.p)
            connect(R_s.n, pos)

            # negative side
            connect(source.p, diode.n)
            connect(source.p, R_sh.n)
            connect(source.p, neg)

            # temperature 
            connect(temp.port, diode.port)
        end
    end

    @mtkbuild sys = Single_Diode()
    return sys
end
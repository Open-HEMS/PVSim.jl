"""
Base model definition for the single diode model.

For a detailed description of the model, see:
https://pvpmc.sandia.gov/modeling-guide/2-dc-module-iv/single-diode-equivalent-circuit-models/

The diode light current is specified using callable parameters:
https://docs.sciml.ai/ModelingToolkit/dev/tutorials/callable_params/#Callable-parameters-and-interpolating-data
"""

"""
    SingleDiode(current, time; n = 1.0, I_s = 1e-6, R_sh_val = 1e5, R_ser_val = 1e-4, T_a = 300.0)

Low level implementation of the single diode model.

TODO: This component should be well-tested.

TODO: Currently the cell temperature is fixed to `T_a`, but this can be made time-dependent.

TODO: Currently the MTK components that make up this model are defined in this file until 
MTK standard libary releases them.

TODO: For now Single_Diode has to be wrapped in a function until some bugs are fixed in MTK.
Most notably, I_L can't be defined inside @components, see: 
https://github.com/SciML/ModelingToolkitStandardLibrary.jl/issues/340

Ultimately we want to set parameters in the model as is shown in: 
https://docs.sciml.ai/ModelingToolkit/dev/basics/MTKLanguage/#@components-begin-block
"""
function SingleDiode(
    current::AbstractVector{<:Real},
    time::AbstractVector{<:Real};
    n::Float64 = 1.0,
    I_s::Float64 = 1e-6,
    R_sh_val::Float64 = 1e5,
    R_ser_val::Float64 = 1e-4,
    T_a::Float64 = 300.0,
)
    @named I_L = ParametrizedInterpolation(LinearInterpolation, current, time)

    @mtkmodel Single_Diode begin

        # TODO: add @parameters block when the bug is fixed

        @components begin
            R_sh = Resistor(R = R_sh_val)
            R_ser = Resistor(R = R_ser_val)
            diode = HeatingDiode(n = n, Is = I_s)
            source = Current()
            temp = FixedTemperature(T = T_a)
            pos = Pin()
            neg = Pin()
        end

        @equations begin
            # current from irradiance
            I_L.input.u ~ t
            connect(I_L.output, source.I)

            # diode anode (+) side connections
            connect(source.n, diode.p)
            connect(source.n, R_sh.p)
            connect(source.n, R_ser.p)

            # diode cathode (-) side connections
            connect(source.p, diode.n)
            connect(source.p, R_sh.n)

            # positive (+) output
            connect(R_ser.n, pos)

            # negative (-) output
            connect(source.p, neg)
            connect(diode.n, neg)
            connect(R_sh.n, neg)

            # temperature 
            connect(temp.port, diode.port)
        end
    end

    @named sys = Single_Diode()
    sys = compose(sys, [I_L])
end


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

"""
    VariableResistor(; name, R_ref = 1.0, T_ref = 300.15, R_const = 1e-3, T_dep = false)

Variable resistor with optional temperature dependency.

The total resistance R ∈ [R_const, R_const + R_ref], where pos is the 
position of the wiper and R_ref is the variable resistance between p and n. 
The total resistance is then:

R = R_const + pos * R_ref

If T_dep is true, then R also depends on the temperature of the heat port with 
temperature coefficient alpha. The total resistance is then:

R = R_const + pos * R_ref * (1 + alpha * (port.T - T_ref))

# States

    - See [OnePort](@ref)
    - `pos(t)`: Position of the wiper (normally 0-1)
    - `R(t)`: Resistance

# Connectors
        
        - `p` Positive pin
        - `n` Negative pin
        - `position` RealInput to set the position of the wiper
        - `port` [HeatPort](@ref) Heat port to model the temperature dependency

# Parameters
    
        - `R_ref`: [`Ω`] Resistance at temperature T_ref when fully closed (pos=1.0)
        - `T_ref`: [K] Reference temperature
        - `R_const`: [`Ω`] Constant resistance between p and n
        - `T_dep`: Temperature dependency
        - `alpha`: [K⁻¹] Temperature coefficient of resistance
        - `enforce_bounds`: Enforce bounds for the position of the wiper (0-1)
"""
@mtkmodel VariableResistor begin
    @extend v, i = oneport = OnePort()

    @structural_parameters begin
        T_dep = false
        enforce_bounds = true
    end

    @parameters begin
        R_ref = 1.0,
        [
            description = "Resistance at temperature T_ref when fully closed (pos=1.0)",
            unit = "Ω",
        ]
        T_ref = 300.15, [description = "Reference temperature", unit = "K"]
        R_const = 1e-3, [description = "Constant resistance between p and n", unit = "Ω"]
    end

    @components begin
        position = RealInput()
    end

    @variables begin
        pos(t), [description = "Position of the wiper (normally 0-1)"]
        R(t), [description = "Resistance", unit = "Ω"]
    end

    if T_dep
        @parameters begin
            alpha =
                1e-3, [description = "Temperature coefficient of resistance", unit = "K^-1"]
        end
        @components begin
            port = HeatPort()
        end
        @equations begin
            port.Q_flow ~ -v * i  # -LossPower
            R ~ R_const + pos * R_ref * (1 + alpha * (port.T - T_ref))
        end
    else
        @equations begin
            R ~ R_const + pos * R_ref
        end
    end

    @equations begin
        pos ~ (enforce_bounds ? clamp(position.u, 0, 1) : position.u)
        v ~ i * R
    end
end
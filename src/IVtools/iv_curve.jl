
function _iv_model(
    light_current::AbstractVector{<:Real},
    time::AbstractVector{<:Real},
    t_end::Float64 = 1.0,
)
    diode = SingleDiode(light_current, time)

    @mtkmodel IV begin
        @parameters begin
            R = 1e5
        end
        @components begin
            load = VariableResistor(R_ref = R, T_dep = false)
            control = Ramp(height = 1.0, duration = t_end)
            ground = Ground()
        end
        @equations begin
            connect(load.position, control.output)
            connect(load.n, ground.g)
            connect(load.p, diode.pos)
            connect(diode.neg, ground.g)
        end
    end

    @named sys = IV()
    sys = compose(sys, [diode])
    sys = structural_simplify(sys)
    return sys
end

function iv_curve(
    light_current::AbstractVector{<:Real},
    time::AbstractVector{<:Real},
    t_end::Float64 = 1.0,
)
    # validate equal length of light_current and time
    if length(light_current) != length(time)
        error(
            "light_current and time must have the same length ($(length(light_current)) != $(length(time)))",
        )
    end

    # create the model
    sys = _iv_model(light_current, time, t_end)
    prob = ODEProblem(sys, [0.0, 0.0], (0.0, t_end))
    sol = solve(prob)
    return sys, sol
end
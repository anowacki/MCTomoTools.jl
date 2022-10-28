# Update a model with a sample from a chain

"""
    update!(m::model, sample::RawSample, sample_index=0) -> model

Using the update specified by `sample`, update `model` so that
the new model is either the same as before if the sample
was rejected, or one of its values is changed according to
the kind of step proposed.

Optionally supply the sample index `sample_index` which is
reported if there are any errors.
"""
function update!(model::Model, sample::RawSample, sample_index=0)
    if Bool(sample.accepted)
        step = Int(sample.step)

        if step === 1
            _update_cell_birth!(model, sample)
        elseif step === 2
            _update_cell_death!(model, sample)
        elseif step === 3
            _update_cell_move!(model, sample)
        elseif step === 4
            _update_velocity!(model, sample)
        elseif step === 5
            _update_body_noise!(model, sample)
        elseif step === 6
            _update_surface_noise!(model, sample)
        elseif step === 7
            _update_event!(model, sample)
        else
            error("unknown step number")
        end
    end
    
    # Check that the model and proposal are consistent
    if Int(sample.ncells) != length(model.nodes)
        _throw_inconsistent_error(sample, model, sample_index)
    end

    nothing
end

function _throw_inconsistent_error(sample, model, sample_index)
    message = "current model (ncells = $(length(model.nodes))) " *
        "and proposal (ncells = $(Int(sample.ncells))) are not consistent" *
        (sample_index == 0 ? "" : " at sample $sample_index")
    error(message)
end

function _update_cell_birth!(m, s)
    push!(m.nodes.coords, @SVector[s.x, s.y, s.z])
    push!(m.nodes.properties, @SVector[s.vp, s.vs, s.density])
    nothing
end

function _update_cell_death!(m, s)
    deleteat!(m.nodes, Int(s.vindex))
    nothing
end

function _update_cell_move!(m, s)
    i = Int(s.vindex)
    m.nodes.coords[i] = @SVector[s.x, s.y, s.z]
    nothing
end

function _update_velocity!(m, s)
    i = Int(s.vindex)
    m.nodes.properties[i] = @SVector[s.vp, s.vs, s.density]
    nothing
end

function _update_body_noise!(m, s)
    i = Int(s.vindex)
    m.body_noise0[i] = s.noise0
    m.body_noise1[i] = s.noise1
    nothing
end

function _update_surface_noise!(m, s)
    i = Int(s.vindex)
    m.surface_noise0[i] = s.noise0
    m.surface_noise1[i] = s.noise1
    nothing
end

function _update_event!(m, s)
    i = Int(s.vindex)
    m.locations[i] = @SVector[s.x, s.y, s.z, s.vp]
    nothing
end

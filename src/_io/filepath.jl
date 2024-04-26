
"""
    get_absolute_filepath(relative_filepath::String)

Returns the absolute filepath by joining the current directory (`@__DIR__`) with the given relative filepath.

# Arguments
- `relative_filepath::String`: The relative filepath to be converted to an absolute filepath.

# Returns
- `String`: The absolute filepath.

"""
function get_absolute_filepath(relative_filepath::String)
    return joinpath(@__DIR__, relative_filepath)
end

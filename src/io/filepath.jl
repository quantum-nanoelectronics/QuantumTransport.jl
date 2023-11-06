function get_absolute_filepath(relative_filepath::String)
    # @__DIR__ gets the directory of the file containing this macro call
    # It is used to ensure that the path is correct no matter where the script is run from
    # joinpath combines @__DIR__ with the relative filepath
    return joinpath(@__DIR__, relative_filepath)
end

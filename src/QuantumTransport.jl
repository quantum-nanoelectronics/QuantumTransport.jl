# The functions here should not be called anywhere else.
# Besides this file, all other items in this directory should be folders containing modules. 

module QuantumTransport


# Given a directory, include all modules in subdirectories
function _helperIncludeModules(dir)
    moduleExists = false
    for item in readdir(dir, join=true)
        if isfile(item) && basename(item) == "Module.jl"
            # println("Including item $item")
            include(item)
            moduleExists = true
        end
    end
    if !moduleExists
        return false
    end
    moduleExistsInAllSubdirs = true
    for item in readdir(dir, join=true)
        if isdir(item) && !_helperIncludeModules(item)
            moduleExistsInAllSubdirs = false
        end
    end
    return moduleExistsInAllSubdirs
end

# Include all modules in subdirectories of the src directory
function _includeModulesInSubdirs(item = @__DIR__)
    # Iterate over all directories in the src directory
    for item in readdir(item, join=true)
        if isdir(item)
            # println("Adding item: $item")
            if !_helperIncludeModules(item)
                error("No module found in directory or subdirectories of $item")
                # println("No module found in directory or subdirectories of $item")
                return false
            end
        end
    end
    return true
end

# Import and export all modules in QuantumTransport
function _importAndExportModules()
    module_names = filter(name -> typeof(getfield(QuantumTransport, name)) <: Module, names(QuantumTransport, all=true))
    for name in module_names
        # Dynamically construct and evaluate the import statement within QuantumTransport
        eval(:(using .$(name)))
        exported_names = names(getfield(QuantumTransport, name), all=false)
        # println("Exported names in $name: ", exported_names)
        for ename in exported_names
            eval(:(export $(ename)))
        end
    end
end

# Call functions in this module
_includeModulesInSubdirs()
_importAndExportModules()

end # module

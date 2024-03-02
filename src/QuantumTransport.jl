# The functions here should not be called anywhere else.
# Besides this file, all other items in this directory should be folders containing modules. 

module QuantumTransport

"""
    _helperIncludeModules(dir)

Function to include modules in all subdirectories of the given directory where all module files are titled Module.jl.

# Arguments
- `dir`: The directory to include modules from.

"""
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


"""
_includeModulesInSubdirs()

This function is responsible for calling _helperIncludeModules(item), where item is each top-level subdirectory of /src. QuantumTransport module.

"""
function _includeModulesInSubdirs()
    # Iterate over all directories in the src directory
    for item in readdir(@__DIR__, join=true)
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
"""
_importAndExportModules()

Dynamically imports and exports modules within the QuantumTransport module.

This function iterates over the modules defined within the QuantumTransport module and dynamically imports them using the `using` statement. It then exports all the names defined within each module using the `export` statement.

No arguments are required for this function.
    
"""
function _importAndExportModules(verbose::Bool=false)
    module_names = filter(name -> typeof(getfield(QuantumTransport, name)) <: Module, names(QuantumTransport, all=true))
    for name in module_names
        # Dynamically construct and evaluate the import statement within QuantumTransport
        eval(:(using .$(name)))
        
        exported_names = names(getfield(QuantumTransport, name), all=false)
        
        # Print information if verbose is true
        if verbose
            println("Exported names in $name: ", join(exported_names, ", "))
        end

        for ename in exported_names
            eval(:(export $(ename)))
        end
    end
end

println("--------LOADING QuantumTransport MODULES----------")

# Call functions in this module
_includeModulesInSubdirs()
_importAndExportModules(true)

println("---------LOADED QuantumTransport MODULES----------")

end # module

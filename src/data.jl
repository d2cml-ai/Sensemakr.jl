
const path = joinpath(dirname(@__FILE__), "..", "data", "darfur.csv");

function load_darfur()

    return CSV.read(path, DataFrame);
end
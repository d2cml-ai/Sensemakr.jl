
const darfur_path = joinpath(dirname(@__FILE__), "..", "data", "darfur.csv");

function load_darfur()

    return CSV.read(darfur_path, DataFrame);
end
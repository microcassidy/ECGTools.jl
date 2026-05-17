using WaveformDB: rdrecord

const DATASET_PATH = ENV["PTB_DATASET"]

function find_headers(root::String =joinpath(DATASET_PATH,"WFDBRecords"))
    map(walkdir(root)) do (path,dirs,files)
        joinpath.(path,filter(x->splitext(x)[2]==".hea",files))
    end |> xs->reduce(vcat,xs)
end

const PTB_HEADERS = find_headers()

export single_signal
function single_signal()
    signal_path = findfirst(endswith("JS00004.hea"),PTB_HEADERS) |> idx -> PTB_HEADERS[idx]
    rdrecord(signal_path,false)
end

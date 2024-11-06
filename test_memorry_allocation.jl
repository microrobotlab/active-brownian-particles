A = rand(1000, 1000)

# Without @views, each slice creates a new array
 for i in 1:10:size(A, 1)
    x=(A[i:i+9, :])
end

# With @views, no new arrays are allocated for slices
@time for i in 1:10:size(A, 1)
    y=(@view A[i:i+9, :])
end

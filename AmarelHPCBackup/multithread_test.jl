# check number of threads
println("Number of threads: $(Threads.nthreads())")

n_arr = [5; 10; 15; 20; 25; 30; 35; 40]

function fibonacci(n::Int64)
    if n <= 1
        return n
    else
        return fibonacci(n - 1) + fibonacci(n - 2)
    end
end

# use multiple threads to compute the fibonacci number for each n in n_arr
fib = zeros(length(n_arr))
Threads.@threads for k = 1:length(n_arr)
    fib[k] = fibonacci(n_arr[k])
end

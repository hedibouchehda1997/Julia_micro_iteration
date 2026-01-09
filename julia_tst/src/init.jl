include("MyLib.jl")
using .MyLib

# Call regular Julia functions (not @ccallable ones)
MyLib.add_numbers(1.0, 2.0)
# MyLib.multiply_arrays([1.0, 2.0], [3.0, 4.0])
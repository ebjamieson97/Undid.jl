using .UNDID
using DataFrames
using StatFiles # For reading .dta


df = DataFrame(load("C:/Users/Eric Bruce Jamieson/Documents/Dalhousie Work/Merit Data Sent from Matt/merit_data_changed.dta"))
columns_to_keep = ["asian", "black", "male", "coll", "merit", "year", "state"]
select!(df, columns_to_keep)
silos = df

names = unique(df.state)
start_times = [minimum(df[df[!, "state"] .== name, "year"]) for name in names]
end_times = [maximum(df[df[!, "state"] .== name, "year"]) for name in names]
treatment_times = [ 
    isempty(df[(df.merit .== 1) .&& (df.state .== name), :year]) ? 0 : minimum(df[(df.merit .== 1) .&& (df.state .== name), :year]) 
    for name in names
]

csv = create_init_csv(names, start_times, end_times, treatment_times)

df = create_diff_df(csv)
# This marks the end of Stage 1 procedures




# Then I need to fill in the empty diff_estimate values at each silo
# (And do some other things)
# Which is Stage 2
df = read_diff_df("C:\\Users\\Eric Bruce Jamieson\\Documents\\Dalhousie Work\\undidjl\\empty_diff_df.csv")


for silo in unique(df[!,"Silo Name"])
    fill_diff_df(silo, df, silos[silos[!, "state"] .== silo,:])
end 



# Stage 3 is then taking all the .csv files produced at each silo and putting them back together 
# Probably need a function to read all the filled .csv files and vcat them together
df = combine_silo_data("C:\\Users\\Eric Bruce Jamieson\\Documents\\Dalhousie Work\\undidjl")

create_ATT_by_gt_by_silo(df, save_as_csv = true)

compute_ATT(create_ATT_by_gt_by_silo(df))


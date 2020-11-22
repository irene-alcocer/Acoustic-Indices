# This script goes over the figure and table captions and creates a proper 
# ordering automatically. This avoids the cumbersome work of having to change
# every figure/table number that comes after a figure/table that was added or 
# removed.
#
# This script assumes that Figure and Tables are labeled as Figure S<number>
# To allow any different naming a proper change must be made in the regex.
#
# First argument is the file path of the Rmd where captions are to be corrected
# Second argument is optional. if we write overwrite as second argument the old
# file will be backed up as <name_file>.backup and the new file will overwrite 
# the old one. Otherwise it will save the new file as <file_name_new>
#
# Example for an index.Rmd file with overwrite option:
#   julia order_captions.jl index.Rmd overwrite

file_path = ARGS[1]
if isassigned(ARGS, 2) 
    overwrite = true
else
    overwrite = false
end

new_file = replace(file_path, r"(.*)(\..*)" => s"\g<1>_new.Rmd")
fw = open(new_file, "w")

lines = open(file_path, "r") do f
    readlines(f, keep = true)
end

fig_nr = 1
table_nr = 1
fig_regex =  r"(.*Figure S)(\d+)(.*)"
table_regex =  r"(.*Table S)(\d+)(.*)"
for line in lines
    if occursin(table_regex, line)
        replaced_number = replace(line, table_regex => 
                                  SubstitutionString("\\g<1>$(table_nr)\\g<3>"))
        global table_nr += 1
        write(fw, replaced_number)
    elseif occursin(fig_regex, line)
        replaced_number = replace(line, fig_regex => 
                                  SubstitutionString("\\g<1>$(fig_nr)\\g<3>"))
        global fig_nr += 1
        write(fw, replaced_number)
    else
        write(fw, line)
    end
end

close(fw)

if overwrite
    mv(file_path, file_path * ".backup")
    mv(new_file, file_path)
end



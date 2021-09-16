# MIT License

# Copyright (c) 2021 Fraunhofer-Chalmers Centre

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.																	

########################################
# Converts the instances from http://soa.iti.es/files/RCmax.7z to csv format used by the LR algorithm
# Download and extract the instances, place this script in the same folder and the new files will be created side-by-side the old files. 
# The instances are also available in the instance.7z, already on the csv format. 
########################################
using DelimitedFiles

cd(@__DIR__ )
folders = [name  for  name in readdir() if isnothing(findfirst('.',name))]
for folder in folders
    for file in readdir(folder)
        A = readdlm(folder*"/"*file)
        B = A[3:end,2:2:end]
        B = collect(transpose(B))
        outfile = replace(file, "txt" => "csv")
        writedlm(folder*"/"*outfile, B, ",")
    end
end
using Dates

struct FuncStruct
    name::String
    func::Function
    # ∇f::Function    # gradient
    # # ∇²::Function{}
    # objFunc::Function
    
end


struct C0Func{t} <: Function
    # name::String
    # func::Function
end

struct C1Func::C0Func
    ∇f::Function
end



# function randIter(func:: Function, exportData::Bool = false; fileName::Union{String,Nothing}=nothing, range::Int=100)
function randIter(func:: Function, exportData::Bool = false; range::Int=100)
    fileName = "ranIter\\"*"randIter_"*Dates.format(Dates.now(),"yyyymmddHHMM")*".txt";    
    println("Func 2");
    if exportData 
        if fileName==nothing
            throw("File name not specified");
        else
            f = open(fileName,"w");
        end
    end
    
    try
        """ Random function to test writing to files"""
        for i ∈ 1:range
            x = func(rand()*i)
            if exportData 
                write(f,"$x\n")
                
            end
        end
        
    finally
        if exportData
            close(f)
        end
    end
    
end 


function randIter(func::Function, file::Union{IOStream, Nothing} = nothing; range::Int = 100)
    """ Random function to test writing to files"""
    println("Func 1")
    for i ∈ 1:range
        x = func(rand()*i)
        if file != nothing && isopen(file)
            write(f,"$x\n")
        end
        
    end
end 



# myFunc(x::Number) = x^2 + 2;
# f = open("test.txt","w");
# fileName = "test.txt";
# randIter(myFunc, true);
# randIter(myFunc, f, range=10);
# close(f);

# Dates.format(Dates.now(),"yyyymmddHHMM");

f(x::Number)::Number = x^2+3;


function newF(x)
    x^2+3;
end


objFunc = FuncStruct("nothing",f);
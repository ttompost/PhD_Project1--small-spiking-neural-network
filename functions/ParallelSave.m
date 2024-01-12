function ParallelSave(SavingName, Val)
% saves a variable from parallel for-loop because save does not work
% author Tea Tompos 21st-Nov-2022
    VarName = inputname(2);
    updateName(VarName, Val);
    save(SavingName,VarName,'-v7.3')
end

function updateName(NewName, Value)
% instead of saving the variable with a name "Val" (which is the name of 2nd input
% variable to ParallelSave), this function captures the true variable name from
% base workspace
% e.g. if variable is called A in base workspace, when passed to
% ParallelSave it will be renamed to Val, but then updateName renames it to
% A again
    assignin('caller',NewName, Value)
end
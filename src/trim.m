%%%This function trims off 1/chopper of the rows from both ends of a
%%%matrix A
function trimmeddata = trim(chopper,A)
A(1:floor(length(A)/chopper),:) = [];
A(floor(length(A)-length(A)/chopper):length(A),:) = [];
trimmeddata = A;
end
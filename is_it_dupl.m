function is_there_a_duplicate = is_it_dupl(A, vect)
%IS_IT_DUPL Returns 1 if column vector vect is already contained as a
%column in matrix A.
is_there_a_duplicate = 0;
for i=1:size(A,2)
    current_col = A(:,i);
    is_it_same = (round(current_col,3) == round(vect,3));
    bool = (is_it_same == ones(size(A,1),1));
    
    if mean(bool) == 1
        is_there_a_duplicate = 1;
        break
    end
end
end


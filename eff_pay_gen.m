clear, clc, tic

% The inputs here are for the Dutch guilder as this currency makes the code run quickly (<10 sec).
% The variables in the workspace that begin with 'an_' are analysis
% variables, which are of particular interest.

% PLEASE NOTE: In line 19 of this code, we have to make a guess about the
% maximum number of unique efficient payments (currently 18). If you get an
% 'out-of-bounds' error, simply increase this number. For currencies with
% multiple different small-valued tokens (like 0.01 and 0.02), we recommend to
% set this number to at least 200.

%% STEP 0: Initialize
D_set = [0.05, 0.1, 0.25, 1, 2.5, 5, 10, 25, 100]'; % denomiantions
A_set = 0.05:0.05:100; % amounts to be taken into consideration
D     = length(D_set);
A     = length(A_set);

P_mat = zeros(D,18,A); 
% P_mat: 3-dimensional array that will hold the efficient payments in its
% columns. Each row corresponds to one token (from D_set, in this order),
% so each cell is a token counter, where negative numbers are change.
% Only the number of non-zero columns will correspond to the number of
% efficient payments for a given amount.
% If you want to analyze this 3-dimensional array, we recommend to loop
% over the third dimensions such that for each amount, we have an efficient
% payment matrix, e.g. for i in A_set: submat = P_mat(:,:,i).

A_set = A_set'; % make it a column vector
round(A_set,3);

%% STEP 1: Cover all elements in A_set that can be paid by only one token
for a = 1:A
    for d = 1:D
        if round(A_set(a),3) == round(D_set(d),3)
            P_mat(d,1,a) = 1;
        end
    end
end
clear d

%% STEP 2: Cover all elements that can be paid with 2 tokens.
for a = 1:A
    if any(P_mat(:,1,a)) == 1 % there is already a single token covering the amount (1st col)
        continue
    else
        n = 1; % counter 
        amount = A_set(a);
        for d1 = 1:D
            for d2 = 1:D
                if round(D_set(d1) + D_set(d2),3) == round(amount,3)
                    P_mat(d1,n,a) = P_mat(d1,n,a) + 1;
                    P_mat(d2,n,a) = P_mat(d2,n,a) + 1; % accounts for case where d1=d2
                    n = n + 1;
                end
                if round(amount - D_set(d1),3) == round(D_set(d2)*(-1),3) % d1: pay, d2: get as change
                    P_mat(d1,n,a) = P_mat(d1,n,a) + 1;
                    P_mat(d2,n,a) = (-1)*P_mat(d2,n,a) - 1; % negative value is change
                    n = n + 1;
                end
            end
        end
        % Remove the duplicates (so duplicate rows) by replacing them by 0.
        submat = P_mat(:,:,a);
        submat = submat'; % 'unique' only works when duplicates are in rows
        [~,ia,~] = unique(submat,'rows','stable');
        i = true(size(submat,1),1);
        i(ia) = false;
        submat(i,:) = 0;
        P_mat(:,:,a) = submat';
    end
end
clear d1 d2 i ia n submat 

prv_tokens = 2; % number of tokens used in the previous step. Here, this is 2.

%% STEP 3
% Idea: Kill one dimension by summing over the absolute values contained in
% all of the respective first columns (do this before and at the end of the while loop). 
% As long as the first column of one 'a' is emty, we have not covered its
% corresposning value yet!
firstcolsums = sum(abs(P_mat(:,1,:)));
temp = permute(firstcolsums,[1 3 2]);
firstcolsums = reshape(temp,[],size(firstcolsums,2),1); % make it a vector
round(firstcolsums,3);

while any(find(firstcolsums==0)) == 1 % as long as there are still first cols with only zeros

    for a = 1:A
        P_submat = P_mat(:,:,a);
        col_sums = sum(abs(P_submat));
        round(col_sums, 3);
        if any(col_sums == prv_tokens) == 0
            continue % skip iteration if there was no efficient number of tokens found in previous step
        end 

        col_indx_set = find(col_sums == round(prv_tokens,3)); % find column indices that belong to a col sum of prv_tokens

        % now, for each col of interest, iteratively add each token with
        % positive and negative sign. Capture this in a Dx2D matrix.
        for j = col_indx_set
            int_col = P_submat(:,j); % column of interest
            add_token_mat = zeros(D,2*D); % matrix to hold the additional tokens

            for k = 1:2*D
                add_token_mat(:,k) = int_col;
            end

            add_token_mat_pos = add_token_mat(:,1:D); % submatrix to hold the positive additions
            add_token_mat_neg = add_token_mat(:,(D+1):(2*D)); % submatrix to hold the negative additions
            for k = 1:D
                add_token_mat_pos(k,k) = add_token_mat_pos(k,k) + 1;
                add_token_mat_neg(k,k) = add_token_mat_neg(k,k) - 1;
            end
            % Now, add_token_mat_neg and add_token_mat_pos hold one additional
            % token in each column. But: Not every column is valid. Only the
            % ones that hold 'prv_tokens+1' tokens as column sums are to be taken
            % into consideration. We thereby account for the case where we 
            % add a token with a positive sign to a combination which has this same
            % token with a negative sign (and vice versa), so this token
            % would cancel out:
            add_token_mat    = [add_token_mat_pos, add_token_mat_neg];
            col_sums_add     = sum(abs(add_token_mat));
            col_indx_set_add = find(col_sums_add == prv_tokens+1);
            add_token_mat    = add_token_mat(:,col_indx_set_add);   
            % Now, add_token_mat only holds columns with a valid amount of
            % tokens. But there is a another restriction: the corresponding
            % amount of each column in add_token_mat needs to be at least as
            % large as the smallest amount in A_set and at most as large as 
            % the largest amount in A_set: 
            col_values       = D_set' * add_token_mat; % case 1: amount too small
            col_indx_set_add = find(col_values >= A_set(1));
            add_token_mat    = add_token_mat(:,col_indx_set_add);

            col_values       = D_set' * add_token_mat; % case 2: amount too big
            col_indx_set_add = find(col_values <= A_set(A));
            add_token_mat    = add_token_mat(:,col_indx_set_add);
            col_values       = round(D_set' * add_token_mat,3);

            % Now: Assign each col in add_token_mat to its corresponding value
            % in P_mat (or P_submat). 

            indx_values = zeros(length(col_values),1);
            for k = 1:length(col_values)
                indx_values(k) = find(round(A_set,3) == round(col_values(k),3));
            end


            for k = 1:length(indx_values) 
                l = indx_values(k); % treat l as the increment            
                P_submat2 = P_mat(:,:,l);
                col_sums2 = sum(abs(P_submat2)); % will contain either 0 or number of tokens
                current_tokens = max(col_sums2); 
                % there are 3 possibilities now: 
                % 1) P_submat2 is empty (if current_tokens == 0) -> append
                % 2) P_submat2 has already been filled with an equally
                % efficient payment (if current_tokens == prvtokens + 1) - > append
                % 3) P_submat2 has already been filled with a more efficient
                % payment (if current_tokens < prvtokens + 1) -> continue

                if round(current_tokens,3) == 0 % case 1)
                    P_submat2(:,1) = add_token_mat(:,k);
                elseif round(current_tokens,3) == round(prv_tokens + 1) % case 2)
                % in this second case, we need to first find out at which col
                % index we can append. Append one col after the largest
                % nonxero col index. But: Append only if this col isn't
                % already contained in the submatrix to avoid duplicates.
                is_it_already_contained = is_it_dupl(P_submat2, add_token_mat(:,k));
                if is_it_already_contained == 1
                    continue % skip iteration if we have this combination already
                end
                
                    nonzero_cols_indx = find(~all(P_submat2==0,1)); % find nonzero columns of matrix
                    larg_nonzero_indx = max(nonzero_cols_indx); % largest nonzero index
                    P_submat2(:,larg_nonzero_indx+1) = add_token_mat(:,k); % append
                else
                    continue
                end % if 
                P_mat(:,:,l) = P_submat2; % submat2 needs to be added back to P_mat
            end % for k
        end % for j
    end % for a
prv_tokens = prv_tokens + 1;

firstcolsums = sum(abs(P_mat(:,1,:)));
temp = permute(firstcolsums,[1 3 2]);
firstcolsums = reshape(temp,[],size(firstcolsums,2),1); % make it a vector


end % while

% P_mat now holds all efficient payments. 

%% ANALYSIS 
tokens_req = zeros(A, prv_tokens);
P_shrunk = P_mat(:,(1:prv_tokens),:);
for a = 1:A
    P_submat = P_shrunk(:,:,a);
    col_sums_add = sum(abs(P_submat));
    tokens_req(a,:) = col_sums_add; % use 'find' to find indices of a given token amount
end
tokens_req_shrunk = tokens_req(:,1);
an_min_tokens = min(tokens_req_shrunk);
an_max_tokens = max(tokens_req_shrunk);
an_med_tokens = median(tokens_req_shrunk);

% for each of the A submatrices, the number of nonzero columns corresponds
% to one efficient payment

an_eff_pay_counter = 0;
an_total_token_counter = 0;
for a = 1:A
    P_submat = P_mat(:,:,a);
    % total number of efficient payments:
    valid_cols = find(~all(P_submat==0)); % cols that are not all zero
    eff_pays_here = length(valid_cols);
    an_eff_pay_counter = an_eff_pay_counter + eff_pays_here;
    % total number of exchanged tokens:
    P_submat = P_submat(:,valid_cols);
    col_sums = sum(abs(P_submat));
    tokens_here = sum(col_sums);
    an_total_token_counter = an_total_token_counter + tokens_here;
end
an_avg_tokens = an_total_token_counter / an_eff_pay_counter ;



clear a add_token_mat add_token_mat_neg add_token_mat_pos col_index__set ...
    col_indx_set_add col_sums col_sums2 col_sums_add col_values current_tokens ...
    eff_pays_here firstcolsums indx_values int_col is_it_already_contained ...
    col_indx_set j k l larg_nonzero_indx P_shrunk P_submat P_submat2 ...
    prv_tokens temp tokens_here nonzero_cols_indx tokens_req tokens_req_shrunk ...
    valid_cols
toc
   
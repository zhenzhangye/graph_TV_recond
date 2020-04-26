function p = backward_solver(K, p, numvex, T, edgeindex,weight)

message_length = 50;
f = K'*T^(-1/2)*p(edgeindex);

K = K./weight;
numedges = size(K,1);
involved_vex = sum(double(sum(abs(K),1)>0));
num_roots = involved_vex - numedges;

involved_vex = find(sum(abs(K),1)>0);

root = zeros(1,num_roots);
root_ind = 0;

visited = ones(1,numvex);
visited(involved_vex) = 0;


tree(numvex,1) = struct('sumA', 0, 'compute_child', 0, 'num_child',0, ...
    'parent',0, 'lambda_pos', 0, 'lambda_neg', 0, ...
    'mes', 0, 'low', 0, 'up', 0, 'children',0);
for i = 1:numvex
    tree(i).mes = zeros(1,message_length);
    tree(i).low = message_length/2;
    tree(i).up = message_length/2;
    tree(i).children = zeros(numvex,1);
    tree(i).num_child = 0;
    tree(i).compute_child = 0;
    tree(i).num_child = 0;
    tree(i).parent = 0;
end

row_used = zeros(1,size(K,1));

while(sum(visited)<numvex)
    root_ind = root_ind+1;
    unvisited_ind = find(visited==0);
    root(root_ind) = unvisited_ind(1);
    queue = zeros(numvex,1);
    front = 1;
    back = 2;
    queue(front) = root(root_ind);
    while(front<back)
        row_id = find(K(:,queue(front))~=0);
        current = queue(front);
        visited(queue(front)) = 1;
        front = front+1;
        for i = 1:length(row_id)
            if(row_used(row_id(i))==0)
            col_id = find(K(row_id(i),:)~=0);
            for j = 1:length(col_id)
                if(visited(col_id(j))==0)
                    queue(back) = col_id(j);
                    back = back+1;
                    tree(col_id(j)).parent = current;
                    tree(current).num_child = tree(current).num_child+1;
                    tree(current).children(tree(current).num_child) = col_id(j);
                end
            end
            row_used(row_id(i))=1;
            end
        end
    end
end

u = zeros(1,numvex);

queue = zeros(numvex,1);
front = 1;
back = 1;
for i = 1:numvex
    if(tree(i).num_child ==0 && sum(double(involved_vex==i))>0)
        queue(back) = i;
        back = back+1;
    end
end

K = K.*weight;
root_computed = 0;

while(front<back)
    current = queue(front);
    front = front +1;
    
    if(~isempty(tree(current).parent) && tree(current).parent>0)
        tree(tree(current).parent).compute_child = tree(tree(current).parent).compute_child+1;
        if(tree(tree(current).parent).compute_child == tree(tree(current).parent).num_child)
            queue(back) = tree(current).parent;
            back = back+1;
        end
    else
        root_computed = root_computed+1;
    end
    
    if(tree(current).num_child~=1)
        tree(current).sumA = 1;
    else
        tree(current).sumA = tree(tree(current).children(1)).sumA +1;
    end
    
    if(tree(current).num_child==0)
        tree(current).low = message_length/2;
        tree(current).up = message_length/2;
        tree(current).mes(message_length/2) = 1-tree(current).sumA;
    else
        if(tree(current).num_child == 1)
            tree(current).low = tree(tree(current).children(1)).low;
            tree(current).up = tree(tree(current).children(1)).up;
            tree(current).mes = tree(tree(current).children(1)).mes;
        else
            tree(current).low = message_length/2;
            tree(current).up = message_length/2;
            tree(current).mes(message_length/2) = 1-tree(current).sumA;
            for i = 1:tree(current).num_child
                mes = zeros(1,message_length);
                child = tree(current).children(i);
                ind = 2;
                indi = tree(current).low;
                indj = tree(child).low;
                mes(1) = tree(current).mes(tree(current).low) + tree(child).mes(tree(child).low) + tree(child).sumA;
                while(indi<tree(current).up && indj<tree(child).up)
                    if(tree(current).mes(indi+1)<tree(child).mes(indj+1))
                        mes(ind) = tree(current).mes(indi+1);
                        ind = ind + 2;
                        mes(ind-1) = tree(current).mes(indi+2)+tree(child).mes(indj)+tree(child).sumA;
                        indi = indi+2;
                    else
                        if(tree(current).mes(indi+1)>tree(child).mes(indj+1))
                            mes(ind) = tree(child).mes(indj+1);
                            ind = ind+2;
                            mes(ind-1) = tree(current).mes(indi) + tree(child).mes(indj+2) + tree(child).sumA;
                            indj = indj+2;
                        else
                            mes(ind) = tree(child).mes(indj+1);
                            mes(ind+1) = tree(child).mes(indj+2)+tree(current).mes(indi+2)+tree(child).sumA;
                            indj = indj+2;
                            ind = ind+2;
                            indi = indi+2;
                        end
                    end
                end
                
                if(indi>=tree(current).up)
                    while(indj<tree(child).up)
                        mes(ind) = tree(child).mes(indj+1);
                        mes(ind+1) = tree(child).mes(indj+2) + tree(current).mes(indi)+tree(child).sumA;
                        indj = indj+2;
                        ind = ind+2;
                    end
                else
                    if(indj>=tree(child).up)
                        while(indi<tree(current).up)
                            mes(ind) = tree(current).mes(indi+1);
                            mes(ind+1) = tree(current).mes(indi+2) + tree(child).mes(indj) + tree(child).sumA;
                            indi = indi+2;
                            ind = ind+2;
                        end
                    end
                end
                
                tree(current).low = message_length/2 - (ind-2)/2;
                tree(current).up = message_length/2 + (ind-2)/2;
                tree(current).mes(tree(current).low:tree(current).up) = mes(1:ind-1);
            end
        end
    end
    
    if(root_computed == num_roots)
        break;
    end
    
    if(length(find(root==current))==0)
    if(tree(current).num_child==0)
        left = -f(current)-weight;
        right = -f(current)+weight;
        tree(current).low = tree(current).low-2;
        tree(current).up = tree(current).up+2;
        tree(current).mes(tree(current).low) = -tree(current).sumA;
        tree(current).mes(tree(current).low+1) = left;
        tree(current).mes(tree(current).up-1) = right;
        tree(current).mes(tree(current).up) = -tree(current).sumA;
        tree(current).lambda_pos = right;
        tree(current).lambda_neg = left;
    else
        child_sum = -weight*tree(current).num_child;
        
        ind = tree(current).low+1;
        
        left = child_sum + tree(current).mes(ind) + f(current);
        
        found = false;
        while(1)
            if(left>-weight)
                found = true;
                break;
            end
            
            if(ind+2>tree(current).up)
                break;
            end
            
            
            
            left = left+(tree(current).sumA + tree(current).mes(ind+1))*(tree(current).mes(ind+2)-tree(current).mes(ind));
            
            ind = ind+2;
            
        end
        
        if(found)
            offset = left-(tree(current).sumA+tree(current).mes(ind-1))*tree(current).mes(ind);
            tree(current).mes(ind-2) = (-weight-offset)/(tree(current).sumA + tree(current).mes(ind-1));
            tree(current).lambda_neg = tree(current).mes(ind-2);
            tree(current).low = ind-3;
            tree(current).mes(ind-3) = -tree(current).sumA;
        else
            offset = left-(tree(current).sumA+tree(current).mes(ind+1))*tree(current).mes(ind);
            tree(current).mes(ind) = (-weight-offset)/(tree(current).sumA + tree(current).mes(ind+1));
            tree(current).lambda_neg = tree(current).mes(ind);
            tree(current).low = ind-1;
            tree(current).mes(ind-1) = -tree(current).sumA;
        end
        
        ind = tree(current).up-1;
        
        if(ind-1==tree(current).low)
            right = -weight;
        else
            right = -child_sum + tree(current).mes(ind) + f(current);
        end
        
        found = false;
        while(1)
            if(right<=weight)
                found = true;
                break;
            end
            
            right = right - (tree(current).sumA+tree(current).mes(ind-1))*(tree(current).mes(ind)-tree(current).mes(ind-2));
            ind = ind-2;
        end
        
        offset = right-(tree(current).sumA + tree(current).mes(ind+1))*tree(current).mes(ind);
        tree(current).lambda_pos = (weight-offset)/(tree(current).sumA + tree(current).mes(ind+1));
        tree(current).mes(ind+2) = tree(current).lambda_pos;
        tree(current).up = ind+3;
        tree(current).mes(ind+3) = -tree(current).sumA;
    end
    end
end

for i = 1:num_roots
    current = root(i);
child_sum = -weight*tree(current).num_child;
ind = tree(current).low+1;
left = child_sum + tree(current).mes(ind)+f(current);
for ind=tree(current).low+1:2:tree(current).up+1
    if(left>0)||(ind>=tree(current).up)
        break;
    end
    left = left + (tree(current).sumA + tree(current).mes(ind+1))*(tree(current).mes(ind+2)-tree(current).mes(ind));
end

offset = left - (tree(current).sumA + tree(current).mes(ind-1))*tree(current).mes(ind);
left = (-offset)/(tree(current).sumA+tree(current).mes(ind-1));
u(current) = left;

queue = zeros(numvex,1);
front = 1;
back = 1;
for i = 1:tree(current).num_child
    queue(back) = tree(current).children(i);
    back = back+1;
end

while(front~=back)
    current = queue(front);
    front = front+1;
    for i = 1:tree(current).num_child
        queue(back) = tree(current).children(i);
        back = back+1;
    end
    
    if(u(tree(current).parent)>=tree(current).lambda_pos)
        u(current) = tree(current).lambda_pos;
    else if(u(tree(current).parent)<=tree(current).lambda_neg)
            u(current) = tree(current).lambda_neg;
        else
            u(current) = u(tree(current).parent);
        end
    end
end
end

p = T^(1/2)*pinv(K')*(u'+f);

end
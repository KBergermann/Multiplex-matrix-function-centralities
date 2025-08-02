%% load string arrays and read out array sizes

names=readtable("names_sci_fi.csv");
layers=table2array(readtable("layers_sci_fi.csv"));

% years 1895, 1904, and 1926 have only one movie with one principal, so no
% edge. Removing these layers...
layers(26,:) = [];  layers(4,:) = []; layers(1,:) = [];

node_id_table = [string(table2array(names(:,1))+ones(size(names(:,1)))), string(table2array(names(:,3)))];

L = size(layers,1);
n = size(node_id_table,1);

%% translate names from "nm0000001"-format to node id's from names_per_title

% bar=waitbar(0,'Looping over layers...');
% (more efficient way?)
parfor (l=1:L, 64)
    opts = detectImportOptions(strcat("title_principals_sci_fi/title_principals_sci_fi_", num2str(layers(l,2)), ".csv"));
    opts.DataLines = [2 Inf];
    opts.VariableTypes(2:end) = {'string'};
    names_table = readtable(strcat("title_principals_sci_fi/title_principals_sci_fi_", num2str(layers(l,2)), ".csv"), opts);
    
    names_per_title = string(table2array(names_table(:,2:end)));
        
    [lines, cols] = size(names_per_title);

    names_per_title_numeric = zeros(lines, cols);
    for i=1:lines
        for j=1:cols
            if ~(ismissing(names_per_title(i,j)))
                try
                    names_per_title_numeric(i,j) = double(node_id_table(find(node_id_table(:,2)==names_per_title(i,j)),1));
                catch
                    if sum(double(node_id_table(:,2)==names_per_title(i,j)))==0
                        fprintf('Warning, skipped node %s in layer %d, which is not in the list of name IDs\n', names_per_title(i,j), l)
                    else
                        % display error message
                        names_per_title_numeric(i,j) = double(node_id_table(find(node_id_table(:,2)==names_per_title(i,j)),1));
                    end
                end
            end
        end
    end

    % build intra-layer adjacency matrix
    % initialize empty block diagonal matrix (adjacency matrix of layer l)
    A_layer = sparse(n,n);

    for i=1:lines
        % insert undirected unweighted edges between all pairs of principals
        nodes = nonzeros(names_per_title_numeric(i,:));
        A_layer(nodes,nodes) = A_layer(nodes,nodes) + 1;
    end

    A_layer = A_layer - diag(diag(A_layer));

    % put layer adjacency on the block diagonal of A_intra
%     if l==1
%         A_intra = A_layer;
%     else
%         A_intra = blkdiag(A_intra, A_layer);
%     end
%     A_intra(n*(l-1)+1:n*l,n*(l-1)+1:n*l) = A_layer;
    A_layer_cell{l} = A_layer;

%     waitbar(l/L,bar);
    fprintf('Completed layer %d (year %d)\n', l, layers(l,2))
end
% close(bar)

A_intra = A_layer_cell{1};
for l=2:L
    A_intra = blkdiag(A_intra, A_layer_cell{l});
end

%% save files
layer_id_table = layers;
save('sci_fi_seasons_intra_layer_adjacency.mat', 'A_intra', 'node_id_table', 'layer_id_table')

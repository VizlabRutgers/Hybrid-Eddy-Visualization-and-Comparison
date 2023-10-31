function [mypaths,mybins,mysizes] = mygetpaths(G)
%
% a path is a connected set of nodes that tracks an eddy (including splits
% and merges) over time

% reorder nodes to hopefully be in sets by path
%[nid,H]=toposort(G);
% hmm, not sure will need this

% this finds connected nodes when multiple sets
% default for DAG requires two directions of connection
% weak type relaxes that
[weak_bins,mysizes] = conncomp(G,'Type','weak');
mybins=unique(weak_bins);

Npaths=max(weak_bins);

% step through nodes adding to path until no connection to previous node
totalnodelist=table2array(G.Nodes);
Nnodes=length(totalnodelist);
%
mypaths=cell(Npaths,1);
for i=1:Nnodes
    curnode=totalnodelist(i);
    nodesbin=weak_bins(i);
    mypaths{nodesbin}=[mypaths{nodesbin}; curnode];
end




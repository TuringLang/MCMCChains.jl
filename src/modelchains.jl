#################### ModelChains ####################

#################### Constructors ####################

function ModelChains(c::Chains, m::Model)
  ModelChains(c.value, c.range, c.names, c.chains, m)
end


#################### Conversions ####################

Base.convert(::Type{Chains}, mc::ModelChains) =
  Chains(mc.value, mc.range, mc.names, mc.chains)


#################### Indexing ####################

function Base.getindex(mc::ModelChains, window, names, chains)
  c = getindex(convert(Chains, mc), window, names2inds(mc, names), chains)
  ModelChains(c, mc.model)
end

names2inds(mc::ModelChains, nodekey::Symbol) = names2inds(mc, [nodekey])

function names2inds(mc::ModelChains, nodekeys::Vector{Symbol})
  inds = Int[]
  missing = Symbol[]
  for key in nodekeys
    keyinds = indexin(names(mc.model, all=key), mc.names)
    0 in keyinds ? push!(missing, key) : append!(inds, keyinds)
  end
  if !isempty(missing)
    throw(ArgumentError(string(
      "chain values are missing for nodes : ",
      join(map(string, missing), ", ")
    )))
  end
  inds
end


function Base.keys(mc::ModelChains, ntype::Symbol, at...)
  values = Symbol[]
  m = mc.model
  nodekeys = ntype == :dependent ?
    keys(m, :dependent) :
    intersect(keys(m, ntype, at...), keys(m, :dependent))
  for key in nodekeys
    all(name -> name in mc.names, names(m, all=key)) && push!(values, key)
  end
  values
end


#################### Auxilliary Functions ####################

function link(c::ModelChains)
  cc = copy(c.value)
  inds_queue = 1:length(c.names)
  for key in intersect(keys(c.model, :monitor), keys(c.model, :stochastic))
    node = c.model[key]
    inds = findall(in(names(node)), c.names)
    if !isempty(inds)
      f(x) = unlist(node, relist(node, x), true)
      cc[:, inds, :] = mapslices(f, cc[:, inds, :], 2)
      inds_queue = setdiff(inds_queue, inds)
    end
  end
  for j in inds_queue
    x = cc[:, j, :]
    if minimum(x) > 0.0
      cc[:, j, :] = maximum(x) < 1.0 ? logit.(x) : log.(x)
    end
  end
  cc
end

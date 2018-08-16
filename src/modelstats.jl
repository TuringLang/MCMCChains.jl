#################### Model-Based Posterior Statistics ####################

function dic(mc::ModelChains)
  nodekeys = keys(mc.model, :output)

  Dhat = -2.0 * logpdf(mc, mean, nodekeys)
  D = -2.0 * logpdf(mc, nodekeys).value
  p = [mean(D) - Dhat, 0.5 * var(D)]

  ChainSummary([Dhat + 2.0 * p  p], ["pD", "pV"],
               ["DIC", "Effective Parameters"], header(mc))
end


function logpdf(mc::ModelChains, f::Function, nodekeys::Vector{Symbol})
  m = mc.model

  relistkeys, updatekeys = getsimkeys(mc, nodekeys)
  relistkeys = union(relistkeys, intersect(nodekeys, keys(m, :block)))
  inds = names2inds(mc, relistkeys)

  m[relistkeys] = relist(m, map(i -> f(mc.value[:, i, :]), inds), relistkeys)
  update!(m, updatekeys)
  mapreduce(key -> logpdf(m[key]), +, nodekeys)
end


logpdf(mc::ModelChains, nodekey::Symbol) = logpdf(mc, [nodekey])

function logpdf(mc::ModelChains,
                nodekeys::Vector{Symbol}=keys(mc.model, :stochastic))
  N = length(mc.range)
  K = size(mc, 3)

  relistkeys, updatekeys = getsimkeys(mc, nodekeys)
  relistkeys = union(relistkeys, intersect(nodekeys, keys(mc.model, :block)))
  inds = names2inds(mc, relistkeys)

  frame = ChainProgressFrame(
    "MCMC Processing of $N Iterations x $K Chain" * "s"^(K > 1), true
  )

  lsts = [
    Any[mc[:, :, k], nodekeys, relistkeys, inds, updatekeys,
        ChainProgress(frame, k, N)]
    for k in 1:K
  ]
  sims = pmap2(logpdf_modelchains_worker, lsts)

  ModelChains(cat(3, sims...), mc.model)
end


function logpdf_modelchains_worker(args::Vector)
  mc, nodekeys, relistkeys, inds, updatekeys, meter = args
  m = mc.model

  sim = Chains(size(mc, 1), 1, start=first(mc), thin=step(mc), names=["logpdf"])

  for i in 1:size(mc.value, 1)
    m[relistkeys] = relist(m, mc.value[i, inds, 1], relistkeys)
    update!(m, updatekeys)
    sim.value[i, 1, 1] = mapreduce(key -> logpdf(m[key]), +, nodekeys)
    next!(meter)
  end

  sim
end


predict(mc::ModelChains, nodekey::Symbol) = predict(mc, [nodekey])

function predict(mc::ModelChains,
                 nodekeys::Vector{Symbol}=keys(mc.model, :output))
  m = mc.model

  outputs = keys(m, :output)
  all(key -> key in outputs, nodekeys) ||
    throw(ArgumentError(string(
      "nodekeys are not all observed Stochastic nodess : ",
      join(map(string, outputs), ", ")
    )))

  nodenames = names(m, nodekeys)
  relistkeys, updatekeys = getsimkeys(mc, nodekeys)
  inds = names2inds(mc, relistkeys)

  c = Chains(size(mc, 1), length(nodenames), chains=size(mc, 3),
             start=first(mc), thin=step(mc), names=nodenames)

  iters, _, chains = size(c.value)
  for k in 1:chains
    for i in 1:iters
      m[relistkeys] = relist(m, mc.value[i, inds, k], relistkeys)
      update!(m, updatekeys)
      f = key -> unlist(m[key], rand(m[key]))
      c.value[i, :, k] = vcat(map(f, nodekeys)...)
    end
  end

  ModelChains(c, m)
end


#################### Auxiliary Functions ####################

function getsimkeys(mc::ModelChains, nodekeys::Vector{Symbol})
  relistkeys = Symbol[]
  updatekeys = Symbol[]

  m = mc.model
  dag = ModelGraph(m)

  nodekeys = intersect(nodekeys, keys(m, :stochastic))
  blockkeys = keys(m, :block)
  dynamickeys = union(blockkeys, keys(m, :target, blockkeys))
  terminalkeys = union(keys(m, :stochastic), keys(mc, :dependent))

  for v in vertices(dag.graph)
    vkey = dag.keys[v]
    if vkey in dynamickeys
      if any(key -> key in nodekeys, gettargets(dag, v, terminalkeys))
        vkey in terminalkeys ?
          push!(relistkeys, vkey) :
          push!(updatekeys, vkey)
      end
    end
  end
  if !isempty(relistkeys) append!(updatekeys, nodekeys) end

  relistkeys, intersect(keys(m, :dependent), updatekeys)
end

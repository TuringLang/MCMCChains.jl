# import Gadfly: draw, Geom, Guide, Layer, layer, PDF, PGF, Plot, plot, PNG, PS,
#        render, Scale, SVG, Theme

#################### Posterior Plots ####################

#################### Generic Methods ####################

# function draw(p::Array{Plot}; fmt::Symbol=:svg, filename::AbstractString="",
#               width::MeasureOrNumber=8inch, height::MeasureOrNumber=8inch,
#               nrow::Integer=3, ncol::Integer=2, byrow::Bool=true,
#               ask::Bool=true)
#
#   fmt in [:pdf, :pgf, :png, :ps, :svg] ||
#     throw(ArgumentError("unsupported draw format $fmt"))
#
#   f(args...) = fmt == :pdf ? PDF(args...) :
#                fmt == :pgf ? PGF(args...) :
#                fmt == :png ? PNG(args...) :
#                fmt == :ps  ? PS(args...)  : SVG(args...)
#
#   isexternalfile = length(filename) > 0
#   addextension = isexternalfile && search(filename, '.') == 0
#   args = (width, height)
#
#   pp = nrow * ncol               ## plots per page
#   ps = length(p)                 ## number of plots
#   np = ceil(Int, ps / pp)        ## number of pages
#
#   mat = Array{Context}(pp)
#   for page in 1:np
#     if ask && page > 1 && !addextension
#       println("Press ENTER to draw next plot")
#       readline(STDIN)
#     end
#
#     if isexternalfile
#       fname = filename
#       if addextension
#         fname = string(fname, '-', page, '.', fmt)
#       end
#       args = (fname, width, height)
#     end
#     img = f(args...)
#
#     nrem = ps - (page - 1) * pp
#     for j in 1:pp
#       if j <= nrem
#         mat[j] = render(p[(page - 1) * pp + j])
#       else
#         mat[j] = context()
#       end
#     end
#     result = byrow ? permutedims(reshape(mat, ncol, nrow), [2, 1]) :
#                      reshape(mat, nrow, ncol)
#
#     draw(img, gridstack(result))
#   end
#
# end

function plot(c::AbstractChains, ptype::Vector{Symbol}=[:trace, :density];
              legend::Bool=false, args...)
  n = length(ptype)
  p = Array{Plot}(n, size(c, 2))
  for i in 1:n
    showlegend = legend && i == n
    p[i, :] = plot(c, ptype[i]; legend=showlegend, args...)
  end
  p
end

function plot(c::AbstractChains, ptype::Symbol; legend::Bool=false, args...)
  ptype == :autocor      ? autocorplot(c; legend=legend, args...) :
  ptype == :bar          ? barplot(c; legend=legend, args...) :
  ptype == :contour      ? contourplot(c; args...) :
  ptype == :density      ? densityplot(c; legend=legend, args...) :
  ptype == :mean         ? meanplot(c; legend=legend, args...) :
  ptype == :mixeddensity ? mixeddensityplot(c; legend=legend, args...) :
  ptype == :trace        ? traceplot(c; legend=legend, args...) :
    throw(ArgumentError("unsupported plot type $ptype"))
end


#################### Plot Engines ####################

function autocorplot(c::AbstractChains;
                     maxlag::Integer=round(Int, 10 * log10(length(c.range))),
                     legend::Bool=false, na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array{Plot}(nvars)
  pos = legend ? :right : :none
  lags = 0:maxlag
  ac = autocor(c, lags=collect(lags))
  for i in 1:nvars
    p = plot(x=repeat(collect(lags * step(c)), outer=[nchains]), vec(ac.value[i, :, :]),
             xaxis="Lag", yaxis="Autocorrelation", title=c.names[i])
    plots[i] = p
  end
  return plots
end

function barplot(c::AbstractChains; legend::Bool=false,
                 position::Symbol=:stack, na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array{Plot}(nvars)
  pos = legend ? :right : :none
  for i in 1:nvars
    # S = unique(c.value[:, i, :])
    # n = length(S)
    # x = repmat(S, 1, nchains)
    # y = zeros(n, nchains)
    # for j in 1:nchains
    #   m = countmap(c.value[:, i, j])
    #   for k in 1:n
    #     if S[k] in keys(m)
    #       y[k, j] = m[S[k]] / nrows
    #     end
    #   end
    # end
    # ymax = maximum(position == :stack ? mapslices(sum, y, 2) : y)
    # plots[i] = plot(x=vec(x), y=vec(y), Geom.bar(position=position),
    #                 color=repeat(c.chains, inner=[n]),
    #                 Scale.color_discrete(), Guide.colorkey(title="Chain"),
    #                 Guide.xlabel("Value", orientation=:horizontal),
    #                 Guide.ylabel("Density", orientation=:vertical),
    #                 Guide.title(c.names[i]), Theme(key_position=pos),
    #                 Scale.y_continuous(minvalue=0.0, maxvalue=ymax))
    p = histogram(vec(c.value[:, i, :]), nbins=50,
                  xaxis="Value", yaxis="Density", title=c.names[i])
    plots[i] = p
  end
  return plots
end

function contourplot(c::AbstractChains; bins::Integer=100, na...)
  nrows, nvars, nchains = size(c.value)
  plots = Plot[]
  offset = 1e4 * eps()
  n = nrows * nchains
  for i in 1:(nvars - 1)
    X = c.value[:, i, :]
    qx = linspace(minimum(X) - offset, maximum(X) + offset, bins + 1)
    mx = map(k -> mean([qx[k], qx[k + 1]]), 1:bins)
    idx = Int[findfirst(k -> qx[k] <= x < qx[k + 1], 1:bins) for x in X]
    for j in (i + 1):nvars
      Y = c.value[:, j, :]
      qy = linspace(minimum(Y) - offset, maximum(Y) + offset, bins + 1)
      my = map(k -> mean([qy[k], qy[k + 1]]), 1:bins)
      idy = Int[findfirst(k -> qy[k] <= y < qy[k + 1], 1:bins) for y in Y]
      density = zeros(bins, bins)
      for k in 1:n
        density[idx[k], idy[k]] += 1.0 / n
      end
      p = contour(mx, my, density, title=c.names[i])
      push!(plots, p)
    end
  end
  return plots
end

function densityplot(c::AbstractChains; legend::Bool=false,
                     trim::Tuple{Real, Real}=(0.025, 0.975), na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array{Plot}(nvars)
  pos = legend ? :right : :none
  for i in 1:nvars
    val = Array{Vector{Float64}}(nchains)
    # for j in 1:nchains
    #   qs = quantile(c.value[:, i, j], [trim[1], trim[2]])
    #   val[j] = c.value[qs[1] .<= c.value[:, i, j] .<= qs[2], i, j]
    # end
    # plots[i] = plot(x=[val...;], Geom.density(),
    #                 color=repeat(c.chains, inner=[length(c.range)]),
    #                 Scale.color_discrete(), Guide.colorkey(title="Chain"),
    #                 Guide.xlabel("Value", orientation=:horizontal),
    #                 Guide.ylabel("Density", orientation=:vertical),
    #                 Guide.title(c.names[i]), Theme(key_position=pos))
    p = histogram(vec(c.value[:, i, :]), nbins=50,
                  xaxis="Value", yaxis="Density", title=c.names[i])
    plots[i] = p
  end
  return plots
end

function meanplot(c::AbstractChains; legend::Bool=false, na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array{Plot}(nvars)
  pos = legend ? :right : :none
  val = cummean(c.value)
  for i in 1:nvars
    p = plot(repeat(collect(c.range), outer=[nchains]), vec(val[:, i, :]),
             xaxis="Iteration", yaxis="Value", title=c.names[i])
    plots[i] = p
  end
  return plots
end

function mixeddensityplot(c::AbstractChains;
                          barbounds::Tuple{Real, Real}=(0, Inf), args...)
  plots = Array{Plot}(size(c, 2))
  discrete = indiscretesupport(c, barbounds)
  plots[discrete] = plot(c[:, discrete, :], :bar; args...)
  plots[.!discrete] = plot(c[:, .!discrete, :], :density; args...)
  return plots
end

function traceplot(c::AbstractChains; legend::Bool=false, na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array{Plot}(nvars)
  pos = legend ? :right : :none
  for i in 1:nvars
    p = plot(repeat(collect(c.range), outer=[nchains]), vec(c.value[:, i, :]),
             xaxis="Iteration", yaxis="Value", title=c.names[i])
    plots[i] = p
  end
  return plots
end

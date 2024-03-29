using Graphs
using MetaGraphs
using GraphMakie
using Makie
using GLMakie
using Makie.Colors
using Makie.ColorSchemes
using NetworkDynamics
using GraphMakie.NetworkLayout: spring
using SciMLBase
using OrderedCollections

export inspect_solution, gparguments

GP_VFAVORITES = [:_ω, :_P, :_Q, :u_arg, :u_mag]
GP_EFAVORITES = [:_P]

function inspect_solution(sol, precord=PRecord(sol.prob))
    # GLMakie.closeall()
    network = _get_nd(sol).f.graph
    fig = Figure(size = (1200, 1200))
    # #####
    # ##### Selectors for sel_nodes and sel_edges
    # #####
    selgrid = fig[1,1] = GridLayout(tellwidth=false, tellheight=true)

    sel_nodes = Observable(OrderedSet{Int}())
    sel_edges = Observable(OrderedSet{Int}())

    selgrid[1,1] = Label(fig, "sel. nodes:")
    selgrid[1,2] = TBSelect(fig, sel_nodes; width=250)
    nplot_btn = selgrid[1,3] = Button(fig, label="node plot", halign=:left)

    selgrid[1,4] = Label(fig, "sel. edges:")
    selgrid[1,5] = TBSelect(fig, sel_edges; width=250)
    eplot_btn = selgrid[1,6] = Button(fig, label="edge plot")

    ####
    #### selector grid to control plot
    ####
    nsymgrid = fig[2,1] = GridLayout(tellwidth=false, tellheight=true)
    nsym_selector = FavSelect(nsymgrid[1,1], GP_VFAVORITES, listvstates(sol, 1:nv(network)); allowmulti=false)
    nstatesym = nsym_selector.selection

    nreltoggle = nsymgrid[1,2] = Toggle(fig)
    nsymgrid[1,3] = Label(fig, "relativ to u0")
    n_rel_to_u0 = nreltoggle.active

    esymgrid = fig[3,1] = GridLayout(tellwidth=false, tellheight=true)
    esym_selector = FavSelect(esymgrid[1,1], GP_EFAVORITES, listestates(sol, 1:ne(network)); allowmulti=false)
    estatesym = esym_selector.selection

    ereltoggle = esymgrid[1,2] = Toggle(fig)
    esymgrid[1,3] = Label(fig, "relativ to u0")
    e_rel_to_u0 = ereltoggle.active

    ####
    #### Statelenses
    ####
    nstatelens = Observable{Any}()
    on(nstatesym; update=true) do sym
        nstatelens[] = vstatef(sol, precord, 1:nv(network), sym; failmode=:warn)
    end
    estatelens = Observable{Any}()
    on(estatesym; update=true) do sym
        estatelens[] = estatef(sol, precord, 1:ne(network), sym; failmode=:warn)
    end

    #####
    ##### Graphplot
    #####
    gpgrid = fig[4,1] = GridLayout()


    ####
    #### Bottom elements
    ####
    bottom_sliders = fig[5,1] = GridLayout(tellwidth=false, tellhight=true)

    ncolorrange, ncolorscheme = _colorrange_slider(bottom_sliders, 1, "node color scaling", nstatelens, n_rel_to_u0)
    ecolorrange, ecolorscheme = _colorrange_slider(bottom_sliders, 2, "edge color scaling", estatelens, e_rel_to_u0)

    Label(bottom_sliders[3,1], text="Time"; halign=:right)
    tslider = Slider(bottom_sliders[3,2], range=Base.range(sol.t[begin], sol.t[end], length=1000))

    t = Observable(0.0)
    connect!(t, tslider.value)
    Label(bottom_sliders[3,3], @lift(@sprintf("%.4f s", $t)))

    tinv_slider = IntervalSlider(bottom_sliders[4,2]; range=tslider.range, tellheight=true)
    Label(bottom_sliders[4,1], @lift(@sprintf("%.4f s", $(tinv_slider.interval)[1])); tellheight=true, halign=:right)
    Label(bottom_sliders[4,3], @lift(@sprintf("%.4f s", $(tinv_slider.interval)[2])); tellheight=true, halign=:left)

    # place the colorbars
    Colorbar(gpgrid[1,1]; colormap=ncolorscheme, colorrange=ncolorrange, vertical=true, label=@lift("node colors "*string($nstatesym[1])), flipaxis=false)
    Colorbar(gpgrid[1,3]; colormap=ecolorscheme, colorrange=ecolorrange, vertical=true, label=@lift("edge colors "*string($estatesym[1])), flipaxis=true)

    ## Graphplot
    gpax = Axis(gpgrid[1,2])

    vstate_vec   = @lift $nstatelens($t)
    vstate_vec_0 = @lift $nstatelens(sol.t[begin])
    node_color = @lift begin
        buf = $vstate_vec
        if ismissing(buf[1])
            buf = zeros(length(buf))
        else
            buf =  buf .- $n_rel_to_u0 * vstate_vec_0[]
        end
        convert(Vector{Float64}, vec(buf))
    end

    estate_vec   = @lift $estatelens($t)
    estate_vec_0 = @lift $estatelens(sol.t[begin])
    edge_color = @lift begin
        buf = $estate_vec
        if ismissing(buf[1])
            buf = zeros(length(buf))
        else
            buf =  buf .- $e_rel_to_u0 * estate_vec_0[]
        end
        convert(Vector{Float64}, vec(buf))
    end

    SMALL = 30
    BIG = 70
    node_size = Observable(fill(SMALL, nv(network)))
    on(sel_nodes; update=true) do selected
        fill!(node_size[], SMALL)
        for sel in selected
            node_size[][sel] = BIG
        end
        notify(node_size)
    end

    THIN = 7
    THICC = 15
    edge_width = Observable(fill(THIN, ne(network)))
    on(sel_edges; update=true) do selected
        fill!(edge_width[], THIN)
        for sel in selected
            edge_width[][sel] = THICC
        end
        notify(edge_width)
    end

    node_marker = let
        markerset = [:circle, :rect, :utriangle, :cross, :diamond, :dtriangle, :pentagon, :xcross]
        groups = sol.prob.f.f.unique_v_indices
        markers = Vector{Symbol}(undef, nv(network))
        for (i,g) in enumerate(groups)
            markers[g] .= markerset[i]
        end
        markers
    end

    graphplot!(gpax,network;
               layout=read_pos_or_spring,
               node_marker,
               node_size,
               node_color,
               node_attr=(;colorrange=ncolorrange, colormap=ncolorscheme),
               edge_width,
               edge_color,
               edge_attr=(;colorrange=ecolorrange, colormap=ecolorscheme));
    hidespines!(gpax)
    hidedecorations!(gpax)

    # delete other interactions on scene
    deregister_interaction!(gpax, :rectanglezoom)

    ####
    #### Click interaction to select and deselct nodes
    ####
    edgeclick = EdgeClickHandler() do idx, event, axis
        sel_edges[] = idx ∈ sel_edges[] ? delete!(sel_edges[], idx) : push!(sel_edges[], idx)
    end
    register_interaction!(gpax, :eclick, edgeclick)

    nodeclick = NodeClickHandler() do idx, event, axis
        sel_nodes[] = idx ∈ sel_nodes[] ? delete!(sel_nodes[], idx) : push!(sel_nodes[], idx)
    end
    register_interaction!(gpax, :nclick, nodeclick)

    ####
    #### Hover interaction to show info
    ####
    HOVER_DEFAULT = "Hover node/edge to see info!"
    hover_text = Observable{String}(HOVER_DEFAULT)
    gpgrid[1, 2] = Label(fig, hover_text, tellwidth=false, tellheight=false, justification=:left, halign=:left, valign=:top, font="JuliaMono")

    nodehover = NodeHoverHandler() do state, idx, event, axis
        if state
            p = NetworkDynamics.p_v_idx(precord(t[]), idx)
            state = vstate_vec[][idx]
            vf = get_vertexf(sol, idx)
            name = vf.name
            d = OrderedDict("State "*string(only(nstatesym[])) => state)
            for (sym, val) in zip(vf.psym, p)
                d[string(sym)] = val
            end
            hover_text[] = "Node $idx: $name\n"*treestyle_string(d)
            node_size[][idx] += 20
            notify(node_size)
        else
            hover_text[] = HOVER_DEFAULT
            notify(sel_nodes)
        end
    end
    register_interaction!(gpax, :nhover, nodehover)

    edgehover = EdgeHoverHandler() do state, idx, event, axis
        if state
            p = NetworkDynamics.p_e_idx(precord(t[]), idx)
            state = estate_vec[][idx]
            # name = get_edgef(sol, idx).name
            d = OrderedDict(first(estatesym[]) => state,
                     "p" => p)
            hover_text[] = "Edge $idx\n"*treestyle_string(d)
            edge_width[][idx] += 4
            notify(edge_width)
        else
            hover_text[] = HOVER_DEFAULT
            notify(sel_edges)
        end
    end
    register_interaction!(gpax, :ehover, edgehover)

    #####
    ##### Keyboard interaction for time
    #####
    register_keyboard_interaction!(fig, tslider)

    #####
    ##### Open node plots in new windows
    #####

    on(nplot_btn.clicks) do n
        f = subplot_window(:node, sol, precord;
                           t=t,
                           tlims=tinv_slider.interval,
                           tslider=tslider,
                           selected=sel_nodes)
        sc = display(GLMakie.Screen(), f)
    end
    on(eplot_btn.clicks) do n
        f = subplot_window(:edge, sol, precord;
                           t=t,
                           tlims=tinv_slider.interval,
                           tslider=tslider,
                           selected=sel_edges)
        sc = display(GLMakie.Screen(), f)
    end

    fig
end

function _colorrange_slider(grid, row, label, statelens, rel)
    Label(grid[row, 1]; justification=:right, text=label)
    sl = Slider(grid[row,2]; range=Base.range(-10, 0, length=100), startvalue=0)
    Label(grid[row,3], @lift(@sprintf("%.2E", 10^$(sl.value))))

    colorscale = @lift 10^$(sl.value)

    ##### Color range stuff
    maxrange = lift(statelens, rel; ignore_equal_values=true) do lens, rel
        values = lens(lens.sol.t)
        if rel
            for col in axes(values,2)
                values[:, 1, col] .-= values[1, 1, col]
            end
        end

        try
            (min, max) = extrema(Iterators.filter(!isnan, skipmissing(values)))
        catch
            (min, max) = (0.0,0.0)
        end
        if min > 0 && max >0
            min = 0.0
        elseif min<0 && max<0
            max = 0.0
        elseif min<0 && max>0
            m = Base.max(-min, max)
            min = -m
            max =  m
        end
        (min, max)
    end

    colorscheme = Observable{ColorScheme}()
    on(maxrange; update=true) do maxrange
        if maxrange == (0.0,0.0)
            colorscheme[] = ColorScheme([colorant"gray80", colorant"gray80", colorant"gray80"])
        elseif (maxrange)[1] < 0
            # ColorScheme([colorant"blue", colorant"gray50", colorant"red"])
            colorscheme[] = ColorSchemes.coolwarm
        else
            colorscheme[] = ColorSchemes.thermal
        end
    end

    colorrange = lift(colorscale, maxrange; ignore_equal_values=true) do colorscale, maxrange
        if (maxrange) == (0.0, 0.0)
            (0.0, 1.0)
        else
            colorscale .* maxrange
        end
    end

    return colorrange, colorscheme
end

function subplot_window(type, sol, precord; t, tlims, tslider, selected)
    fig = Figure(size=(1000, 800))
    symgrid = fig[1,1] = GridLayout(tellwidth=false, tellheight=true)

    sym_selector = if type === :node
        FavSelect(symgrid[1,1], GP_VFAVORITES, @lift(listvstates(sol, $selected)); allowmulti=true)
    elseif type === :edge
        FavSelect(symgrid[1,1], GP_EFAVORITES, @lift(listestates(sol, $selected)); allowmulti=true)
    end
    statesyms = sym_selector.selection

    reltoggle = symgrid[1,2] = Toggle(fig)
    symgrid[1,3] = Label(fig, "relativ to u0")
    rel_to_u0 = reltoggle.active

    ax = Axis(fig[2, 1])
    on(tlims, update=true) do lims
        xlims!(ax, lims)
    end

    statelens = Observable{Any}()
    if type === :node
        onany(statesyms, selected) do syms, idxs
            # @info "update node statelens"
            statelens[] = vstatef(sol, precord, idxs, syms; failmode=:warn)
        end
    elseif type === :edge
        onany(statesyms, selected) do syms, idxs
            # @info "update edge statelens"
            statelens[] = estatef(sol, precord, idxs, syms; failmode=:warn)
        end
    end

    ts = Observable(range(sol.t[begin], sol.t[end], length=1000))
    on(statelens) do _
        ts.val = range(sol.t[begin], sol.t[end], length=1000)
    end

    lastupdate = Ref(time())
    on(ax.finallimits) do lims
        lastupdate[] = time()
    end

    timer = Timer(.5; interval=.5) do _
        if time() > lastupdate[] + 0.5
            lims = ax.finallimits[]
            tmin = max(sol.t[begin], lims.origin[1])
            tmax = min(sol.t[end], tmin + lims.widths[1])
            ts.val = range(tmin, tmax, length=1000)
            lastupdate[] = Inf
            notify(rel_to_u0)
        end
    end

    on(events(fig.scene).window_open) do open
        if !open
            @info "Stop Timer"
            close(timer)
        end
    end


    data = Observable{Array{Union{Missing, Float32}}}(zeros(1,1,1))
    onany(statelens, rel_to_u0) do lens, rel
        @info "Resample"
        newdat = lens(ts[])
        if rel
            for idxidx in axes(newdat , 3), symidx in axes(newdat , 2)
                newdat[:, symidx, idxidx] .-= newdat[begin, symidx, idxidx]
            end
        end
        # @info "after" data[] newdat
        if size(data[]) !== size(newdat)
            data.val = newdat
            empty!(data.listeners)
            # TODO: maybe delete listeners for data here?
        else
            data[] = newdat
        end
    end

    legendref = Ref{Legend}()
    on(statelens; priority=-1) do _
        @info "Replot"
        # if isassigned(legendref)
        #     delete!(legendref[])
        # end
        empty!(ax)
        vlines!(ax, t; color=:black)
        for idxidx in axes(data[], 3)
            for symidx in axes(data[], 2)
                # TODO: this is bad because series stays around
                series = @lift view($data, :, symidx, idxidx)
                ismissing(series[][begin]) && continue
                lines!(ax, ts, series;
                       label=string(statesyms[][symidx])*(rel_to_u0[] ? " (rel)" : "")* " @ "*string(selected[][idxidx]),
                       color=Cycled(idxidx),
                       linestyle=ax.scene.theme.palette.linestyle[][symidx])
            end
        end
        # if !isempty(data[])
        #    legendref[] = axislegend(ax)
        # end
    end

    ####
    #### Click interaction to set time
    ####
    set_time_interaction = (event::MouseEvent, axis) -> begin
        if event.type === MouseEventTypes.leftclick
            pos = mouseposition(axis.scene)[1]
            set_close_to!(tslider, pos)
            return Consume(true)
        end
        return Consume(false)
    end
    register_interaction!(set_time_interaction, ax, :set_time)

    # arrows to move time
    register_keyboard_interaction!(fig, tslider)

    notify(selected) # shoud trigger everything
    return fig
end

function register_keyboard_interaction!(fig, tslider)
    on(fig.scene.events.keyboardbutton) do e
        if e.action == Keyboard.press || e.action == Keyboard.repeat
            if e.key == Keyboard.left
                steps = shift_pressed(fig) ? 10 : 1
                mv_slider(tslider, -steps)
                return true
            elseif e.key == Keyboard.right
                steps = shift_pressed(fig) ? 10 : 1
                mv_slider(tslider, steps)
                return true
            end
        end
        return false
    end
end

function shift_pressed(fig)
    keys = fig.scene.events.keyboardstate
    Makie.Keyboard.left_shift ∈ keys || Makie.Keyboard.right_shift ∈ keys
end

function mv_slider(slider, steps=1)
    idx = slider.selected_index[]
    set_close_to!(slider, slider.range[][clamp(idx+steps, firstindex(slider.range[]), lastindex(slider.range[]))])
end

function gparguments(sol::ODESolution,
                     precord::PRecord,
                     network;
                     t,
                     nstatelens,
                     ncolorrange,
                     ncolorscheme,
                     n_rel_to_u0,
                     sel_nodes = Observable(Set{Int}()))

    t::Observable = t isa Observable ? t : Observable(t)
    NV = nv(network)

    node_marker = try
        markdict = Dict(:load => :rect, :gen => :circle, :syncon => :circle);
        [markdict[k] for k in get_prop(network, 1:NV, :type)];
    catch
        markerset = [:circle, :rect, :utriangle, :cross, :diamond, :dtriangle, :pentagon, :xcross]
        groups = sol.prob.f.f.unique_v_indices
        markers = Vector{Symbol}(undef, nv(network))
        for (i,g) in enumerate(groups)
            markers[g] .= markerset[i]
        end
        markers
    end

    u0statevec = Observable(Vector{Float64}(undef, NV))

    on(nstatelens; update=true) do lens
        u0statevec[] .= @views lens(sol.t)[1, :]
    end

    statevec = Observable(Vector{Float64}(undef, NV))

    onany(t, nstatelens, n_rel_to_u0) do t, lens, rel
        statevec[] .= @views lens(t)[1, :]
        if rel
            statevec[] .-= u0statevec[]
        end
        notify(statevec)
    end

    node_color = Observable(Vector{RGB{Float64}}(undef, NV))
    onany(statevec, ncolorrange, ncolorscheme) do statevec, range, scheme
        for i in 1:NV
            node_color[][i] = isnan(statevec[i]) ? RGB(0,0,0) : get(scheme, statevec[i], range)
        end
        notify(node_color)
    end
    notify(t)

    SMALL = 30
    BIG = 50
    node_size = Observable(fill(SMALL, NV))
    onany(sel_nodes) do selected
        fill!(node_size[], SMALL)
        for sel in selected
            node_size[][sel] = BIG
        end
        notify(node_size)
    end
    notify(sel_nodes)

    return (;layout=read_pos_or_spring,
            node_marker,
            node_color,
            node_size,
            node_attr=(;colorrange=ncolorrange, colormap=ncolorscheme));
end

read_pos_or_spring(g) = spring(g)
function read_pos_or_spring(g::MetaGraph)
    pos = get_prop(g, 1:nv(g), :pos)
    if any(ismissing.(pos))
        return spring(g)
    else
        return pos
    end
end

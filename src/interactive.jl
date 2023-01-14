using Graphs
using Printf
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

function inspect_solution(sol, network=sol.prob.f.f.graph, precord=PRecord(sol.prob))
    fig = Figure(resolution = (1200, 1200))
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
        hover_text[] = if state
            p = NetworkDynamics.p_v_idx(precord(t[]), idx)
            state = vstate_vec[][idx]
            d = OrderedDict(first(nstatesym[]) => state,
                     "p" => p)
            "Node $idx\n"*treestyle_string(d)
        else
            HOVER_DEFAULT
        end
    end
    register_interaction!(gpax, :nhover, nodehover)

    edgehover = EdgeHoverHandler() do state, idx, event, axis
        hover_text[] = if state
            p = NetworkDynamics.p_e_idx(precord(t[]), idx)
            state = estate_vec[][idx]
            d = OrderedDict(first(nstatesym[]) => state,
                     "p" => p)
            "Edge $idx\n"*treestyle_string(d)
        else
            HOVER_DEFAULT
        end
    end
    register_interaction!(gpax, :ehover, edgehover)

    #####
    ##### Keyboard interaction for time
    #####
    register_keyboard_interaction!(fig, tslider)
    return fig

    #####
    ##### Open node plots in new windows
    #####
    on(nselectors[3].clicks) do n
        f = nodeplot_window(sol, precord, tslider, sel_nodes; tlims=tinterval_slider.interval)
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
            (min, max) = extrema(skipmissing(values))
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

function subplot_window(type, sol, precord, tslider, tlims, selected)
    fig = Figure(resolution=(1000, 800))

    if type === :node
        symgrid = fig[2,1] = GridLayout(tellwidth=false, tellheight=true)
        sym_selector = FavSelect(symgrid[1,1], GP_VFAVORITES, listvstates(sol, 1:nv(network)); allowmulti=false)
        nstatesym = sym_selector.selection

        nreltoggle = symgrid[1,2] = Toggle(fig)
        symgrid[1,3] = Label(fig, "relativ to u0")
        n_rel_to_u0 = nreltoggle.active
    elseif type === :edge
        esymgrid = fig[3,1] = GridLayout(tellwidth=false, tellheight=true)
        esym_selector = FavSelect(esymgrid[1,1], GP_EFAVORITES, listestates(sol, 1:ne(network)); allowmulti=false)
        estatesym = esym_selector.selection

        ereltoggle = esymgrid[1,2] = Toggle(fig)
        esymgrid[1,3] = Label(fig, "relativ to u0")
        e_rel_to_u0 = ereltoggle.active
    end


    ax = Axis(fig[3, 1])

    plots = Dict{Int, Any}()

    onany(syms, n_rel_to_u0) do syms, rel
        @debug "Syms chagned to $syms "*(rel ? "(relative) " : "") * "clear all."
        empty!(ax)
        empty!(plots)
        Makie.vlines!(ax, tslider.value; color=:black)
    end
    notify(syms)

    on(tlims, update=true) do lims
        xlims!(ax, lims)
    end

    legend = nothing
    onany(sel_nodes, syms, n_rel_to_u0) do selected, syms, rel
        added   = setdiff(selected, keys(plots))
        removed = setdiff(keys(plots), selected)
        for i in added
            plist = []
            for (isym, s) in enumerate(syms)
                try
                    ts = timeseries(sol, precord, i, s)
                    if rel
                        ts = TimeSeries(ts.t, ts.x .- ts.x[begin], ts.name)
                    end
                    p = lines!(ax, ts;
                               label=string(s)*(rel ? " (rel)" : "")* " @ "*string(i),
                               linewidth=3,
                               color=Cycled(i),
                               linestyle=ax.palette.linestyle[][isym])
                    push!(plist, p)
                    @debug "Added plot $s "*(rel ? " (rel)" : "")*" for node $i"
                catch e
                    @debug "Could not plot $s "*(rel ? " (rel)" : "")*" for node $i" e
                    # variable not found
                end
            end
            plots[i] = plist
        end
        for i in removed
            plist = plots[i]
            delete!(plots, i)
            for p in plist
                delete!(ax, p)
            end
        end
        #= TODO: Legend broken
        if !(isempty(added) && isempty(removed))
            if !isnothing(legend)
                # legend.width[] = 0
                # legend.height[] = 0
                # legend.tellwidth[]=false
                # delete!(legend)
            end
            # isnothing(legend) || delete!(legend)
            if !isempty(plots)
                # legend = Legend(fig[2,2], ax, ":"*string(syms)*" at nodes:"; bgcolor=:white)
            end
            # display(GLMakie.Screen(fig.scene), fig)
        end
        =#
    end

    # initialize box
    symbox.displayed_string[] = string(syms[][1])
    symbox.stored_string[] = string(syms[][1])

    # arrows to move time
    register_keyboard_interaction!(fig, tslider)

    # click to set time
    # TODO: this blocks the rectagle zoom interaction
    set_time_interaction = (event::MouseEvent, axis) -> begin
        if event.type === MouseEventTypes.rightclick
            pos = mouseposition(axis.scene)[1]
            set_close_to!(tslider, pos)
            return true
        end
        return false
    end
    register_interaction!(set_time_interaction, ax, :set_time)

    fig
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

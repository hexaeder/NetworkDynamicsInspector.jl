using Makie

Makie.@Block TBSelect begin
    @forwarded_layout
    textbox::Textbox
    selection::Observable{<:Any}
    @attributes begin
        "Controls if multiple syms can be selected at the same time."
        allowmulti = true
        "The horizontal alignment of the block in its suggested bounding box."
        halign = :center
        "The vertical alignment of the block in its suggested bounding box."
        valign = :center
        "The width setting of the block."
        width = Auto()
        "The height setting of the block."
        height = Auto()
        "Controls if the parent layout can adjust to this block's width"
        tellwidth::Bool = false
        "Controls if the parent layout can adjust to this block's height"
        tellheight::Bool = true
        "The align mode of the block in its parent GridLayout."
        alignmode = Inside()
    end
end
function Makie.initialize_block!(ts::TBSelect, selection=Observable(Symbol[]))
    setfield!(ts, :selection, selection)
    CT = eltype(selection)
    ET = eltype(CT)

    ts.textbox = Textbox(ts.layout[1,1], width=ts.width)

    ts.textbox.stored_string.ignore_equal_values = true
    on(ts.textbox.stored_string) do str
        str = replace(str, r"[\s+]"=>s"") # remove whitespace
        if occursin(r"^\s*$", str) # empty
            ts.selection[] = CT()
        elseif !ts.allowmulti[]
            ts.selection[] = CT([_convert_text(ET, str)])
        else
            parts = split(str, ','; keepempty=false)
            ts.selection[] = CT(_convert_text.(ET, parts))
        end
    end

    on(ts.selection; update=true) do sel
        rep = repr_selection(sel)
        try
            ts.textbox.displayed_string[] = isempty(rep) ? " " : rep
        catch e
            # try again, weird makie bug
            ts.textbox.displayed_string[] = isempty(rep) ? " " : rep
        end
    end
    # ts.textbox.stored_string.val = rep
end
_convert_text(ET::Type{Any}, str) = str
_convert_text(ET::Type{String}, str) = str
_convert_text(ET::Type{Symbol}, str) = Symbol(replace(str,r"^:"=>s""))
_convert_text(ET::Type{<:Number}, str) = parse(ET, str)
function repr_selection(sel)
    isempty(sel) && return ""
    str = reduce(*, string.(sel) .* ", ")
    first(str, length(str)-2)
end


Makie.@Block FavSelect begin
    @forwarded_layout
    buttons::Vector{Button}
    menu::Menu
    textbox::TBSelect
    selection::Observable{Vector{Symbol}}
    @attributes begin
        "Controls if multiple syms can be selected at the same time."
        allowmulti = true
        "The horizontal alignment of the block in its suggested bounding box."
        halign = :center
        "The vertical alignment of the block in its suggested bounding box."
        valign = :center
        "The width setting of the block."
        width = Auto()
        "The height setting of the block."
        height = Auto()
        "Controls if the parent layout can adjust to this block's width"
        tellwidth::Bool = false
        "Controls if the parent layout can adjust to this block's height"
        tellheight::Bool = true
        "The align mode of the block in its parent GridLayout."
        alignmode = Inside()
    end
end

function Makie.initialize_block!(fs::FavSelect, favs, syms)
    syms = syms isa Observable ? syms : Observable(syms)
    setfield!(fs, :selection, Observable([favs[1]]))

    fs.buttons = Button[]
    for (i, fav) in enumerate(favs)
        btn = Button(fs.layout[1,i], label=string(fav))
        push!(fs.buttons, btn)
    end
    fs.menu = Menu(fs.layout[1,length(favs)+1], options=syms)

    fs.textbox = TBSelect(fs.layout[1,length(favs)+2], fs.selection; width=150, allowmulti=fs.allowmulti)

    for i in 1:length(favs)
        on(fs.buttons[i].clicks) do n
            _update_selection(fs, favs[i])
        end
    end

    on(fs.menu.selection) do sel
        isnothing(sel) || _update_selection(fs, sel)
    end

    on(fs.selection) do sel
        # handle menu
        sel_from_menu = intersect(sel, syms[])
        if isempty(sel_from_menu)
            fs.menu.prompt[] = "nothing selected"
        else
            fs.menu.prompt[] = repr_selection(sel_from_menu) * " selected"
        end
        fs.menu.i_selected[] = 0

        # buttons
        for (i, fav) in enumerate(favs)
            if fav ∈ sel
                fs.buttons[i].buttoncolor = Makie.COLOR_ACCENT[]
            else
                fs.buttons[i].buttoncolor = RGBf(0.94,0.94,0.94)
            end
        end
    end

    notify(fs.selection)

    return nothing
end

function _update_selection(fs::FavSelect, sym)
    if fs.allowmulti[] && shift_pressed(fs)
        if sym ∉ fs.selection[]
            push!(fs.selection[], sym)
            notify(fs.selection)
        end
    else
        fs.selection[] = [sym]
    end
end

shift_pressed(x::Makie.Block) = shift_pressed(x.parent)
shift_pressed(x::Figure) = shift_pressed(x.scene)
function shift_pressed(sc::Scene)
    keys = sc.events.keyboardstate
    Makie.Keyboard.left_shift ∈ keys || Makie.Keyboard.right_shift ∈ keys
end

